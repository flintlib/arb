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
_fmprb_poly_rgamma_series(fmprb_ptr res, fmprb_ptr h, long hlen, long len, long prec)
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

        gamma_lgamma_series_at_one(u, len, wp);

        _fmprb_vec_neg(u, u, len);
        _fmprb_poly_exp_series(t, u, len, len, wp);

        if (r == 1)
        {
            _fmprb_vec_swap(v, t, len);
        }
        else if (r <= 0)
        {
            fmprb_set(f, h);
            fmprb_one(f + 1);
            rflen = FLINT_MIN(len, 2 - r);
            _fmprb_poly_rfac_series_ui(u, f, FLINT_MIN(2, len), 1 - r, rflen, wp);
            _fmprb_poly_mullow(v, t, len, u, rflen, len, wp);
        }
        else
        {
            fmprb_one(f);
            fmprb_one(f + 1);
            rflen = FLINT_MIN(len, r);
            _fmprb_poly_rfac_series_ui(v, f, FLINT_MIN(2, len), r - 1, rflen, wp);

            /* TODO: use div_series? */
            _fmprb_poly_inv_series(u, v, rflen, len, wp);
            _fmprb_poly_mullow(v, t, len, u, len, len, wp);
        }
    }
    else
    {
        /* otherwise use Stirling series */
        gamma_stirling_choose_param_fmprb(&reflect, &r, &n, h, 1, 0, wp);

        /* rgamma(h) = (gamma(1-h+r) sin(pi h)) / (rf(1-h, r) * pi), h = h0 + t*/
        if (reflect)
        {
            /* u = gamma(r+1-h) */
            fmprb_sub_ui(f, h, r + 1, wp);
            fmprb_neg(f, f);
            gamma_stirling_eval_fmprb_series(t, f, n, len, wp);
            _fmprb_poly_exp_series(u, t, len, len, wp);
            for (i = 1; i < len; i += 2)
                fmprb_neg(u + i, u + i);

            /* v = sin(pi x) */
            fmprb_const_pi(f + 1, wp);
            fmprb_mul(f, h, f + 1, wp);
            _fmprb_poly_sin_series(v, f, 2, len, wp);

            _fmprb_poly_mullow(t, u, len, v, len, len, wp);

            /* rf(1-h,r) * pi */
            if (r == 0)
            {
                fmprb_const_pi(u, wp);
                _fmprb_vec_scalar_div(v, t, len, u, wp);
            }
            else
            {
                fmprb_sub_ui(f, h, 1, wp);
                fmprb_neg(f, f);
                fmprb_set_si(f + 1, -1);
                rflen = FLINT_MIN(len, r + 1);
                _fmprb_poly_rfac_series_ui(v, f, FLINT_MIN(2, len), r, rflen, wp);
                fmprb_const_pi(u, wp);
                _fmprb_vec_scalar_mul(v, v, rflen, u, wp);

                /* divide by rising factorial */
                /* TODO: might better to use div_series, when it has a good basecase */
                _fmprb_poly_inv_series(u, v, rflen, len, wp);
                _fmprb_poly_mullow(v, t, len, u, len, len, wp);
            }
        }
        else
        {
            /* rgamma(h) = rgamma(h+r) rf(h,r) */
            if (r == 0)
            {
                fmprb_add_ui(f, h, r, wp);
                gamma_stirling_eval_fmprb_series(t, f, n, len, wp);
                _fmprb_vec_neg(t, t, len);
                _fmprb_poly_exp_series(v, t, len, len, wp);
            }
            else
            {
                fmprb_set(f, h);
                fmprb_one(f + 1);
                rflen = FLINT_MIN(len, r + 1);
                _fmprb_poly_rfac_series_ui(t, f, FLINT_MIN(2, len), r, rflen, wp);

                fmprb_add_ui(f, h, r, wp);
                gamma_stirling_eval_fmprb_series(v, f, n, len, wp);
                _fmprb_vec_neg(v, v, len);
                _fmprb_poly_exp_series(u, v, len, len, wp);

                _fmprb_poly_mullow(v, u, len, t, rflen, len, wp);
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
fmprb_poly_rgamma_series(fmprb_poly_t res, const fmprb_poly_t f, long n, long prec)
{
    if (f->length == 0 || n == 0)
    {
        fmprb_poly_zero(res);
    }
    else
    {
        fmprb_poly_fit_length(res, n);
        _fmprb_poly_rgamma_series(res->coeffs, f->coeffs, f->length, n, prec);
        _fmprb_poly_set_length(res, n);
        _fmprb_poly_normalise(res);
    }
}

