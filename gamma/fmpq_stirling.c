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

#include "gamma.h"

void
gamma_fmpq_stirling(fmprb_t y, const fmpq_t a, long prec)
{
    int reflect;
    long r, n, wp;
    fmprb_t x, t, u, v;

    wp = prec + FLINT_BIT_COUNT(prec);

    fmprb_init(x);
    fmprb_init(t);
    fmprb_init(u);
    fmprb_init(v);

    fmprb_set_fmpq(x, a, wp);
    gamma_stirling_choose_param_fmprb(&reflect, &r, &n, x, 1, 0, wp);

    if (reflect)
    {
        /* gamma(x) = (rf(1-x, r) * pi) / (gamma(1-x+r) sin(pi x)) */
        fmpq_t b;
        fmpq_init(b);
        fmpz_sub(fmpq_numref(b), fmpq_denref(a), fmpq_numref(a));
        fmpz_set(fmpq_denref(b), fmpq_denref(a));
        gamma_rising_fmprb_fmpq_ui_bsplit(u, b, r, wp);
        fmpq_clear(b);
        fmprb_const_pi(v, wp);
        fmprb_mul(u, u, v, wp);
        fmprb_sub_ui(t, x, 1, wp);
        fmprb_neg(t, t);
        fmprb_add_ui(t, t, r, wp);
        gamma_stirling_eval_series_fmprb(v, t, n, 0, wp);
        fmprb_exp(v, v, wp);
        fmprb_sin_pi_fmpq(t, a, wp);
        fmprb_mul(v, v, t, wp);
    }
    else
    {
        /* gamma(x) = gamma(x+r) / rf(x,r) */
        fmprb_add_ui(t, x, r, wp);
        gamma_stirling_eval_series_fmprb(u, t, n, 0, wp);
        fmprb_exp(u, u, prec);
        gamma_rising_fmprb_fmpq_ui_bsplit(v, a, r, wp);
    }

    fmprb_div(y, u, v, prec);

    fmprb_clear(t);
    fmprb_clear(u);
    fmprb_clear(v);
    fmprb_clear(x);
}

