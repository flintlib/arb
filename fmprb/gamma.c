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

    Copyright (C) 2012, 2013 Fredrik Johansson

******************************************************************************/

#include "fmprb.h"
#include "gamma.h"

static void
_fmprb_gamma(fmprb_t y, const fmprb_t x, long prec, int inverse)
{
    int reflect;
    long r, n, wp;
    fmprb_t t, u, v;

    if (fmprb_is_exact(x))
    {
        const fmpr_struct * mid = fmprb_midref(x);

        if (fmpr_is_special(mid))
        {
            if (!inverse && fmpr_is_pos_inf(mid))
            {
                fmprb_set(y, x);
            }
            else if (fmpr_is_nan(mid) || fmpr_is_neg_inf(mid) || !inverse)
            {
                fmpr_nan(fmprb_midref(y));
                fmpr_pos_inf(fmprb_radref(y));
            }
            else
            {
                fmprb_zero(y);
            }
            return;
        }
        else
        {
            /* fast gamma(n), gamma(n/2) or gamma(n/4) */
            if (fmpr_cmpabs_2exp_si(mid, prec) < 0 &&
                fmpr_is_int_2exp_si(mid, -2))
            {
                fmpq_t a;
                fmpq_init(a);
                fmpr_get_fmpq(a, mid);
                fmprb_gamma_fmpq(y, a, prec + 2 * inverse);
                if (inverse)
                    fmprb_inv(y, y, prec);
                fmpq_clear(a);
                return;
            }
        }
    }

    wp = prec + FLINT_BIT_COUNT(prec);

    gamma_stirling_choose_param_fmprb(&reflect, &r, &n, x, 1, 0, wp);

    fmprb_init(t);
    fmprb_init(u);
    fmprb_init(v);

    if (reflect)
    {
        /* gamma(x) = (rf(1-x, r) * pi) / (gamma(1-x+r) sin(pi x)) */
        fmprb_sub_ui(t, x, 1, wp);
        fmprb_neg(t, t);
        gamma_rising_fmprb_ui_bsplit(u, t, r, wp);
        fmprb_const_pi(v, wp);
        fmprb_mul(u, u, v, wp);
        fmprb_add_ui(t, t, r, wp);
        gamma_stirling_eval_fmprb(v, t, n, 0, wp);
        fmprb_exp(v, v, wp);
        fmprb_sin_pi(t, x, wp);
        fmprb_mul(v, v, t, wp);
    }
    else
    {
        /* gamma(x) = gamma(x+r) / rf(x,r) */
        fmprb_add_ui(t, x, r, wp);
        gamma_stirling_eval_fmprb(u, t, n, 0, wp);
        fmprb_exp(u, u, prec);
        gamma_rising_fmprb_ui_bsplit(v, x, r, wp);
    }

    if (inverse)
        fmprb_div(y, v, u, prec);
    else
        fmprb_div(y, u, v, prec);

    fmprb_clear(t);
    fmprb_clear(u);
    fmprb_clear(v);
}

void
fmprb_gamma(fmprb_t y, const fmprb_t x, long prec)
{
    _fmprb_gamma(y, x, prec, 0);
}

void
fmprb_rgamma(fmprb_t y, const fmprb_t x, long prec)
{
    _fmprb_gamma(y, x, prec, 1);
}

void
fmprb_lgamma(fmprb_t y, const fmprb_t x, long prec)
{
    int reflect;
    long r, n, wp;
    fmprb_t t, u;

    wp = prec + FLINT_BIT_COUNT(prec);

    gamma_stirling_choose_param_fmprb(&reflect, &r, &n, x, 0, 0, wp);

    /* log(gamma(x)) = log(gamma(x+r)) - log(rf(x,r)) */
    fmprb_init(t);
    fmprb_init(u);

    fmprb_add_ui(t, x, r, wp);
    gamma_stirling_eval_fmprb(u, t, n, 0, wp);
    gamma_rising_fmprb_ui_bsplit(t, x, r, wp);
    fmprb_log(t, t, wp);
    fmprb_sub(y, u, t, prec);

    fmprb_clear(t);
    fmprb_clear(u);
}

