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

void
_fmprb_poly_pow_series(fmprb_ptr h,
    fmprb_srcptr f, long flen,
    fmprb_srcptr g, long glen, long len, long prec)
{
    if (glen == 1)
    {
        _fmprb_poly_pow_fmprb_series(h, f, flen, g, len, prec);
        return;
    }

    /* f^g = exp(g * log(f)) */
    if (flen == 1)
    {
        fmprb_t t;
        fmprb_init(t);
        fmprb_log(t, f, prec);
        _fmprb_vec_scalar_mul(h, g, glen, t, prec);
        _fmprb_poly_exp_series(h, h, glen, len, prec);
        fmprb_clear(t);
    }
    else
    {
        fmprb_ptr t;
        t = _fmprb_vec_init(len);
        _fmprb_poly_log_series(t, f, flen, len, prec);
        _fmprb_poly_mullow(h, t, len, g, glen, len, prec);
        _fmprb_poly_exp_series(h, h, len, len, prec);
        _fmprb_vec_clear(t, len);

    }
}

void
fmprb_poly_pow_series(fmprb_poly_t h,
    const fmprb_poly_t f, const fmprb_poly_t g, long len, long prec)
{
    long flen, glen;

    flen = f->length;
    glen = g->length;

    flen = FLINT_MIN(flen, len);
    glen = FLINT_MIN(glen, len);

    if (len == 0)
    {
        fmprb_poly_zero(h);
        return;
    }

    if (glen == 0)
    {
        fmprb_poly_one(h);
        return;
    }

    if (flen == 0)
    {
        fmprb_poly_zero(h);
        return;
    }

    if (flen == 1 && glen == 1)
    {
        fmprb_poly_fit_length(h, 1);
        fmprb_pow(h->coeffs, f->coeffs, g->coeffs, prec);
        _fmprb_poly_set_length(h, 1);
        _fmprb_poly_normalise(h);
        return;
    }

    if (f == h || g == h)
    {
        fmprb_poly_t t;
        fmprb_poly_init2(t, len);
        _fmprb_poly_pow_series(t->coeffs,
            f->coeffs, flen, g->coeffs, glen, len, prec);
        _fmprb_poly_set_length(t, len);
        _fmprb_poly_normalise(t);
        fmprb_poly_swap(t, h);
        fmprb_poly_clear(t);
    }
    else
    {
        fmprb_poly_fit_length(h, len);
        _fmprb_poly_pow_series(h->coeffs,
            f->coeffs, flen, g->coeffs, glen, len, prec);
        _fmprb_poly_set_length(h, len);
        _fmprb_poly_normalise(h);
    }
}

