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
#include "zeta.h"

void
_fmprb_poly_gamma_series(fmprb_ptr res, fmprb_ptr h, long hlen, long len, long prec)
{
    int reflect;
    long i, rflen, r, n, wp;
    fmprb_ptr t, u, v;
    fmprb_struct f[2];

    hlen = FLINT_MIN(hlen, len);
    wp = prec + FLINT_BIT_COUNT(prec);

    t = _fmprb_vec_init(len);
    u = _fmprb_vec_init(len);
    v = _fmprb_vec_init(len);
    fmprb_init(f);
    fmprb_init(f + 1);

    /* use zeta values at small integers */
    if (fmprb_is_int(h) && (fmpr_cmpabs_ui(fmprb_midref(h), prec / 2) < 0))
    {
        r = fmpr_get_si(fmprb_midref(h), FMPR_RND_DOWN);

        if (r <= 0)
        {
            _fmprb_vec_indeterminate(v, len);
        }
        else
        {
            fmprb_zero(u);
            if (len > 1) fmprb_const_euler(u + 1, wp);
            if (len > 2) zeta_ui_vec(u + 2, 2, len - 2, wp);
            for (i = 2; i < len; i++)
                fmprb_div_ui(u + i, u + i, i, wp);
            for (i = 1; i < len; i += 2)
                fmprb_neg(u + i, u + i);

            if (r == 1)
            {
                _fmprb_poly_exp_series(v, u, len, len, wp);
            }
            else
            {
                _fmprb_poly_exp_series(t, u, len, len, wp);
                fmprb_one(f);
                fmprb_one(f + 1);
                rflen = FLINT_MIN(len, r);
                _fmprb_poly_rfac_series_ui(u, f, FLINT_MIN(2, len), r - 1, rflen, wp);
                _fmprb_poly_mullow(v, t, len, u, rflen, len, wp);
            }
        }
    }
    else
    {
        /* otherwise use Stirling series */
        gamma_stirling_choose_param_fmprb(&reflect, &r, &n, h, 1, 0, wp);

        /* gamma(h) = (rf(1-h, r) * pi) / (gamma(1-h+r) sin(pi h)), h = h0 + t*/
        if (reflect)
        {
            /* u = 1/gamma(r+1-h) */
            fmprb_sub_ui(f, h, r + 1, wp);
            fmprb_neg(f, f);
            gamma_stirling_eval_fmprb_series(t, f, n, len, wp);
            _fmprb_vec_neg(t, t, len);
            _fmprb_poly_exp_series(u, t, len, len, wp);
            for (i = 1; i < len; i += 2)
                fmprb_neg(u + i, u + i);

            /* v = 1/sin(pi x) */
            fmprb_const_pi(f + 1, wp);
            fmprb_mul(f, h, f + 1, wp);
            _fmprb_poly_sin_series(t, f, 2, len, wp);
            _fmprb_poly_inv_series(v, t, len, len, wp);

            _fmprb_poly_mullow(t, u, len, v, len, len, wp);

            /* rf(1-h,r) * pi */
            if (r == 0)
            {
                rflen = 1;
                fmprb_const_pi(u, wp);
            }
            else
            {
                fmprb_sub_ui(f, h, 1, wp);
                fmprb_neg(f, f);
                fmprb_set_si(f + 1, -1);
                rflen = FLINT_MIN(len, r + 1);
                _fmprb_poly_rfac_series_ui(u, f, FLINT_MIN(2, len), r, rflen, wp);
                fmprb_const_pi(v, wp);
                _fmprb_vec_scalar_mul(u, u, rflen, v, wp);
            }

            /* multiply by rising factorial */
            _fmprb_poly_mullow(v, t, len, u, rflen, len, wp);
        }
        else
        {
            /* gamma(h) = gamma(h+r) / rf(h,r) */
            if (r == 0)
            {
                fmprb_add_ui(f, h, r, wp);
                gamma_stirling_eval_fmprb_series(t, f, n, len, wp);
                _fmprb_poly_exp_series(v, t, len, len, wp);
            }
            else
            {
                /* TODO: div_series may be better (once it has a good basecase),
                         if the rising factorial is short */
                fmprb_set(f, h);
                fmprb_one(f + 1);
                rflen = FLINT_MIN(len, r + 1);
                _fmprb_poly_rfac_series_ui(u, f, FLINT_MIN(2, len), r, rflen, wp);
                _fmprb_poly_inv_series(t, u, rflen, len, wp);

                fmprb_add_ui(f, h, r, wp);
                gamma_stirling_eval_fmprb_series(v, f, n, len, wp);
                _fmprb_poly_exp_series(u, v, len, len, wp);

                _fmprb_poly_mullow(v, u, len, t, len, len, wp);
            }
        }
    }

    /* compose with nonconstant part */
    fmprb_zero(t);
    _fmprb_vec_set(t + 1, h + 1, hlen - 1);
    _fmprb_poly_compose_series(res, v, len, t, hlen, len, prec);

    fmprb_clear(f);
    fmprb_clear(f + 1);
    _fmprb_vec_clear(t, len);
    _fmprb_vec_clear(u, len);
    _fmprb_vec_clear(v, len);
}

void
fmprb_poly_gamma_series(fmprb_poly_t res, const fmprb_poly_t f, long n, long prec)
{
    fmprb_poly_fit_length(res, n);

    if (f->length == 0 || n == 0)
        _fmprb_vec_indeterminate(res->coeffs, n);
    else
        _fmprb_poly_gamma_series(res->coeffs, f->coeffs, f->length, n, prec);

    _fmprb_poly_set_length(res, n);
    _fmprb_poly_normalise(res);
}

