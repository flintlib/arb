/*=============================================================================

    This file is part of ARB.

    ARB is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    ARB is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ARB; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include "fmprb_poly.h"

static __inline__ void
_fmprb_poly_get_scale(fmpz_t scale, fmprb_srcptr x, long xlen, fmprb_srcptr y, long ylen)
{
    long xa, xb, ya, yb, den;

    fmpz_zero(scale);

    /* ignore zeros (and infs/nans!); find the first and last
       finite nonzero entries to determine the scale */
    xa = 0;
    xb = xlen - 1;
    while (xa < xlen && fmpr_is_special(fmprb_midref(x + xa))) xa++;
    while (xb > xa && fmpr_is_special(fmprb_midref(x + xb))) xb--;

    ya = 0;
    yb = ylen - 1;
    while (ya < ylen && fmpr_is_special(fmprb_midref(y + ya))) ya++;
    while (yb > ya && fmpr_is_special(fmprb_midref(y + yb))) yb--;

    /* compute average of exponent differences, weighted by the lengths */
    if (xa <= xb && ya <= yb && (xa < xb || ya < yb))
    {
        fmpz_add(scale, scale, fmpr_expref(fmprb_midref(x + xb)));
        fmpz_add_ui(scale, scale, fmpz_bits(fmpr_manref(fmprb_midref(x + xb))));
        fmpz_sub(scale, scale, fmpr_expref(fmprb_midref(x + xa)));
        fmpz_sub_ui(scale, scale, fmpz_bits(fmpr_manref(fmprb_midref(x + xa))));

        fmpz_add(scale, scale, fmpr_expref(fmprb_midref(y + yb)));
        fmpz_add_ui(scale, scale, fmpz_bits(fmpr_manref(fmprb_midref(y + yb))));
        fmpz_sub(scale, scale, fmpr_expref(fmprb_midref(y + ya)));
        fmpz_sub_ui(scale, scale, fmpz_bits(fmpr_manref(fmprb_midref(y + ya))));

        den = (xb - xa) + (yb - ya);

        /* scale = floor(scale / den + 1/2) = floor((2 scale + den) / (2 den)) */
        fmpz_mul_2exp(scale, scale, 1);
        fmpz_add_ui(scale, scale, den);
        fmpz_fdiv_q_ui(scale, scale, 2 * den);
    }
}

/* copy exponents deeply, mantissas shallowly */
static __inline__ void
_fmpr_set_shallow_mul2exp(fmpr_t y, const fmpr_t x, const fmpz_t e)
{
    *fmpr_manref(y) = *fmpr_manref(x);

    if (fmpr_is_special(x))
    {
        fmpz_init_set(fmpr_expref(y), fmpr_expref(x));
    }
    else
    {
        fmpz_init(fmpr_expref(y));
        fmpz_add_inline(fmpr_expref(y), fmpr_expref(x), e);
    }
}

static __inline__ void
_fmprb_set_shallow_mul2exp(fmprb_t y, const fmprb_t x, const fmpz_t e)
{
    _fmpr_set_shallow_mul2exp(fmprb_midref(y), fmprb_midref(x), e);
    _fmpr_set_shallow_mul2exp(fmprb_radref(y), fmprb_radref(x), e);
}

static __inline__ void
_fmprb_clear_shallow(fmprb_t x)
{
    fmpz_clear(fmpr_expref(fmprb_midref(x)));
    fmpz_clear(fmpr_expref(fmprb_radref(x)));
}

void
_fmprb_poly_mullow_block_scaled(fmprb_ptr z,
    fmprb_srcptr x, long xlen,
    fmprb_srcptr y, long ylen,
    long len, long prec)
{
    fmpz_t scale;
    fmpz_init(scale);
    _fmprb_poly_get_scale(scale, x, xlen, y, ylen);

    if (fmpz_is_zero(scale))
    {
        _fmprb_poly_mullow_block(z, x, xlen, y, ylen, len, prec);
    }
    else
    {
        fmprb_ptr xtmp, ytmp;
        fmpz_t t;
        long i;

        xtmp = flint_malloc(sizeof(fmprb_struct) * (xlen + ylen));
        ytmp = xtmp + xlen;

        fmpz_init(t);
        fmpz_zero(t);
        for (i = 0; i < FLINT_MAX(xlen, ylen); i++)
        {
            if (i < xlen)
                _fmprb_set_shallow_mul2exp(xtmp + i, x + i, t);
            if (i < ylen)
                _fmprb_set_shallow_mul2exp(ytmp + i, y + i, t);
            fmpz_sub(t, t, scale);
        }

        _fmprb_poly_mullow_block(z, xtmp, xlen, ytmp, ylen, len, prec);

        fmpz_zero(t);
        for (i = 0; i < len; i++)
        {
            fmprb_mul_2exp_fmpz(z + i, z + i, t);
            fmpz_add(t, t, scale);
        }

        for (i = 0; i < xlen; i++) _fmprb_clear_shallow(xtmp + i);
        for (i = 0; i < ylen; i++) _fmprb_clear_shallow(ytmp + i);

        fmpz_clear(t);
        flint_free(xtmp);
    }

    fmpz_clear(scale);
}

void
fmprb_poly_mullow_block_scaled(fmprb_poly_t res, const fmprb_poly_t poly1,
              const fmprb_poly_t poly2, long n, long prec)
{
    long xlen, ylen, zlen;

    xlen = poly1->length;
    ylen = poly2->length;

    if (xlen == 0 || ylen == 0 || n == 0)
    {
        fmprb_poly_zero(res);
        return;
    }

    xlen = FLINT_MIN(xlen, n);
    ylen = FLINT_MIN(ylen, n);
    zlen = FLINT_MIN(xlen + ylen - 1, n);

    if (res == poly1 || res == poly2)
    {
        fmprb_poly_t tmp;
        fmprb_poly_init2(tmp, zlen);
        _fmprb_poly_mullow_block_scaled(tmp->coeffs, poly1->coeffs, xlen,
            poly2->coeffs, ylen, zlen, prec);
        fmprb_poly_swap(res, tmp);
        fmprb_poly_clear(tmp);
    }
    else
    {
        fmprb_poly_fit_length(res, zlen);
        _fmprb_poly_mullow_block_scaled(res->coeffs, poly1->coeffs, xlen,
            poly2->coeffs, ylen, zlen, prec);
    }

    _fmprb_poly_set_length(res, zlen);
    _fmprb_poly_normalise(res);
}

