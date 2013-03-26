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

#include "fmpcb.h"
#include "gamma.h"

void
fmpcb_digamma(fmpcb_t y, const fmpcb_t x, long prec)
{
    int reflect;
    long r, n, wp;
    fmpcb_t t, u, v;

    wp = prec + FLINT_BIT_COUNT(prec);

    gamma_stirling_choose_param_fmpcb(&reflect, &r, &n, x, 1, 1, wp);

    fmpcb_init(t);
    fmpcb_init(u);
    fmpcb_init(v);

    /* psi(x) = psi((1-x)+r) - h(1-x,r) - pi*cot(pi*x) */
    if (reflect)
    {
        fmpcb_sub_ui(t, x, 1, wp);
        fmpcb_neg(t, t);
        fmpcb_cot_pi(v, x, wp);
        fmprb_const_pi(fmpcb_realref(u), wp);
        fmpcb_mul_fmprb(v, v, fmpcb_realref(u), wp);
        gamma_harmonic_sum_fmpcb_ui_bsplit(u, t, r, wp);
        fmpcb_add(v, v, u, wp);
        fmpcb_add_ui(t, t, r, wp);
        gamma_stirling_eval_series_fmpcb(u, t, n, 1, wp);
        fmpcb_sub(y, u, v, wp);
    }
    else
    {
        fmpcb_add_ui(t, x, r, wp);
        gamma_stirling_eval_series_fmpcb(u, t, n, 1, wp);
        gamma_harmonic_sum_fmpcb_ui_bsplit(t, x, r, wp);
        fmpcb_sub(y, u, t, prec);
    }

    fmpcb_clear(t);
    fmpcb_clear(u);
    fmpcb_clear(v);
}

