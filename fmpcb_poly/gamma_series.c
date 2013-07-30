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

#include "fmpcb_poly.h"
#include "gamma.h"
#include "zeta.h"

void
_fmpcb_poly_gamma_series(fmpcb_ptr res, fmpcb_srcptr h, long hlen, long len, long prec)
{
    int reflect;
    long i, rflen, r, n, wp;
    fmpcb_ptr t, u, v;
    fmpcb_struct f[2];

    hlen = FLINT_MIN(hlen, len);
    wp = prec + FLINT_BIT_COUNT(prec);

    t = _fmpcb_vec_init(len);
    u = _fmpcb_vec_init(len);
    v = _fmpcb_vec_init(len);
    fmpcb_init(f);
    fmpcb_init(f + 1);

    /* TODO: use real code at real numbers */
    if (0)
    {
    }
    else
    {
        /* otherwise use Stirling series */
        gamma_stirling_choose_param_fmpcb(&reflect, &r, &n, h, 1, 0, wp);

        /* gamma(h) = (rf(1-h, r) * pi) / (gamma(1-h+r) sin(pi h)), h = h0 + t*/
        if (reflect)
        {
            /* u = 1/gamma(r+1-h) */
            fmpcb_sub_ui(f, h, r + 1, wp);
            fmpcb_neg(f, f);
            gamma_stirling_eval_fmpcb_series(t, f, n, len, wp);
            _fmpcb_vec_neg(t, t, len);
            _fmpcb_poly_exp_series(u, t, len, len, wp);
            for (i = 1; i < len; i += 2)
                fmpcb_neg(u + i, u + i);

            /* v = 1/sin(pi x) */
            fmpcb_const_pi(f + 1, wp);
            fmpcb_mul(f, h, f + 1, wp);
            _fmpcb_poly_sin_series(t, f, 2, len, wp);
            _fmpcb_poly_inv_series(v, t, len, len, wp);

            _fmpcb_poly_mullow(t, u, len, v, len, len, wp);

            /* rf(1-h,r) * pi */
            if (r == 0)
            {
                rflen = 1;
                fmpcb_const_pi(u, wp);
            }
            else
            {
                fmpcb_sub_ui(f, h, 1, wp);
                fmpcb_neg(f, f);
                fmpcb_set_si(f + 1, -1);
                rflen = FLINT_MIN(len, r + 1);
                _fmpcb_poly_rising_ui_series(u, f, FLINT_MIN(2, len), r, rflen, wp);
                fmpcb_const_pi(v, wp);
                _fmpcb_vec_scalar_mul(u, u, rflen, v, wp);
            }

            /* multiply by rising factorial */
            _fmpcb_poly_mullow(v, t, len, u, rflen, len, wp);
        }
        else
        {
            /* gamma(h) = gamma(h+r) / rf(h,r) */
            if (r == 0)
            {
                fmpcb_add_ui(f, h, r, wp);
                gamma_stirling_eval_fmpcb_series(t, f, n, len, wp);
                _fmpcb_poly_exp_series(v, t, len, len, wp);
            }
            else
            {
                /* TODO: div_series may be better (once it has a good basecase),
                         if the rising factorial is short */
                fmpcb_set(f, h);
                fmpcb_one(f + 1);
                rflen = FLINT_MIN(len, r + 1);
                _fmpcb_poly_rising_ui_series(u, f, FLINT_MIN(2, len), r, rflen, wp);
                _fmpcb_poly_inv_series(t, u, rflen, len, wp);

                fmpcb_add_ui(f, h, r, wp);
                gamma_stirling_eval_fmpcb_series(v, f, n, len, wp);
                _fmpcb_poly_exp_series(u, v, len, len, wp);

                _fmpcb_poly_mullow(v, u, len, t, len, len, wp);
            }
        }
    }

    /* compose with nonconstant part */
    fmpcb_zero(t);
    _fmpcb_vec_set(t + 1, h + 1, hlen - 1);
    _fmpcb_poly_compose_series(res, v, len, t, hlen, len, prec);

    fmpcb_clear(f);
    fmpcb_clear(f + 1);
    _fmpcb_vec_clear(t, len);
    _fmpcb_vec_clear(u, len);
    _fmpcb_vec_clear(v, len);
}

void
fmpcb_poly_gamma_series(fmpcb_poly_t res, const fmpcb_poly_t f, long n, long prec)
{
    fmpcb_poly_fit_length(res, n);

    if (f->length == 0 || n == 0)
        _fmpcb_vec_indeterminate(res->coeffs, n);
    else
        _fmpcb_poly_gamma_series(res->coeffs, f->coeffs, f->length, n, prec);

    _fmpcb_poly_set_length(res, n);
    _fmpcb_poly_normalise(res);
}

