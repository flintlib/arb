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
#include "gamma.h"

void
_fmprb_poly_lgamma_series(fmprb_ptr res, fmprb_ptr h, long hlen, long len, long prec)
{
    int reflect;
    long r, n, rflen, wp;
    fmprb_t zr;
    fmprb_ptr t, u, f;

    hlen = FLINT_MIN(hlen, len);
    wp = prec + FLINT_BIT_COUNT(prec);
    gamma_stirling_choose_param_fmprb(&reflect, &r, &n, h, 0, 0, wp);

    t = _fmprb_vec_init(len);
    u = _fmprb_vec_init(len);
    f = _fmprb_vec_init(2);
    fmprb_init(zr);

    /* log(gamma(z)) = log(gamma(z+r)) - log(rf(z,r)) */
    fmprb_add_ui(zr, h, r, wp);

    /* u = log(gamma(z+r)) */
    gamma_stirling_eval_fmprb_series(u, zr, n, len, wp);

    /* subtract t = log(rf(z,r)) */
    if (r != 0)
    {
        fmprb_set(f, h);
        fmprb_one(f + 1);
        rflen = FLINT_MIN(len, r + 1);
        _fmprb_poly_rfac_series_ui(t, f, FLINT_MIN(2, len), r, rflen, prec);
        _fmprb_poly_log_series(t, t, rflen, len, wp);
        _fmprb_vec_sub(u, u, t, len, prec);
    }

    /* compose with nonconstant part */
    fmprb_zero(t);
    _fmprb_vec_set(t + 1, h + 1, hlen - 1);
    _fmprb_poly_compose_series(res, u, len, t, hlen, len, prec);

    fmprb_clear(zr);
    _fmprb_vec_clear(t, len);
    _fmprb_vec_clear(u, len);
    _fmprb_vec_clear(f, 2);
}

void
fmprb_poly_lgamma_series(fmprb_poly_t res, const fmprb_poly_t f, long n, long prec)
{
    if (f->length == 0 || n == 0)
    {
        printf("fmprb_poly_lgamma_series: require n > 0 and nonzero input\n");
        abort();
    }

    fmprb_poly_fit_length(res, n);
    _fmprb_poly_lgamma_series(res->coeffs, f->coeffs, f->length, n, prec);
    _fmprb_poly_set_length(res, n);
    _fmprb_poly_normalise(res);
}

