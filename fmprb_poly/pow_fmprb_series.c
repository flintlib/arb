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
_fmprb_poly_pow_fmprb_series(fmprb_ptr h,
    fmprb_srcptr f, long flen, const fmprb_t g, long len, long prec)
{
    int f_binomial, g_exact, g_int;

    while (flen > 0 && fmprb_is_zero(f + flen - 1))
        flen--;

    if (flen <= 1)
    {
        fmprb_pow(h, f, g, prec);
        _fmprb_vec_zero(h + 1, len - 1);
        return;
    }

    g_exact = fmprb_is_exact(g);
    g_int = fmprb_is_int(g);
    f_binomial = _fmprb_vec_is_zero(f + 1, flen - 2);

    /* g = small integer */
    if (g_exact && g_int &&
            fmpr_cmpabs_2exp_si(fmprb_midref(g), FLINT_BITS - 1) < 0)
    {
        long e, hlen;

        e = fmpz_get_si(fmpr_manref(fmprb_midref(g))) <<
            fmpz_get_ui(fmpr_expref(fmprb_midref(g)));

        hlen = poly_pow_length(flen, FLINT_ABS(e), len);

        if (e >= 0)
        {
            _fmprb_poly_pow_ui_trunc_binexp(h, f, flen, e, hlen, prec);
            _fmprb_vec_zero(h + hlen, len - hlen);
            return;
        }
        else if (!f_binomial)
        {
            fmprb_ptr t;
            t = _fmprb_vec_init(hlen);
            _fmprb_poly_pow_ui_trunc_binexp(t, f, flen, -e, hlen, prec);
            _fmprb_poly_inv_series(h, t, hlen, len, prec);
            _fmprb_vec_clear(t, hlen);
            return;
        }
    }

    /* (a + bx^c)^g */
    if (f_binomial)
    {
        long i, j, d;
        fmprb_t t;

        fmprb_init(t);

        d = flen - 1;
        fmprb_pow(h, f, g, prec);
        fmprb_div(t, f + d, f, prec);

        for (i = 1, j = d; j < len; i++, j += d)
        {
            fmprb_sub_ui(h + j, g, i - 1, prec);
            fmprb_mul(h + j, h + j, h + j - d, prec);
            fmprb_mul(h + j, h + j, t, prec);
            fmprb_div_ui(h + j, h + j, i, prec);
        }

        if (d > 1)
        {
            for (i = 1; i < len; i++)
                if (i % d != 0)
                    fmprb_zero(h + i);
        }

        fmprb_clear(t);
        return;
    }

    /* g = +/- 1/2 */
    if (g_exact && fmpr_cmpabs_2exp_si(fmprb_midref(g), -1) == 0)
    {
        if (fmpr_sgn(fmprb_midref(g)) > 0)
            _fmprb_poly_sqrt_series(h, f, flen, len, prec);
        else
            _fmprb_poly_rsqrt_series(h, f, flen, len, prec);
        return;
    }

    /* f^g = exp(g*log(f)) */
    _fmprb_poly_log_series(h, f, flen, len, prec);
    _fmprb_vec_scalar_mul(h, h, len, g, prec);
    _fmprb_poly_exp_series(h, h, len, len, prec);

}

void
fmprb_poly_pow_fmprb_series(fmprb_poly_t h,
    const fmprb_poly_t f, const fmprb_t g, long len, long prec)
{
    long flen;

    flen = f->length;
    flen = FLINT_MIN(flen, len);

    if (len == 0)
    {
        fmprb_poly_zero(h);
        return;
    }

    if (fmprb_is_zero(g))
    {
        fmprb_poly_one(h);
        return;
    }

    if (flen == 0)
    {
        fmprb_poly_zero(h);
        return;
    }

    if (f == h)
    {
        fmprb_poly_t t;
        fmprb_poly_init2(t, len);
        _fmprb_poly_pow_fmprb_series(t->coeffs, f->coeffs, flen, g, len, prec);
        _fmprb_poly_set_length(t, len);
        _fmprb_poly_normalise(t);
        fmprb_poly_swap(t, h);
        fmprb_poly_clear(t);
    }
    else
    {
        fmprb_poly_fit_length(h, len);
        _fmprb_poly_pow_fmprb_series(h->coeffs, f->coeffs, flen, g, len, prec);
        _fmprb_poly_set_length(h, len);
        _fmprb_poly_normalise(h);
    }
}

