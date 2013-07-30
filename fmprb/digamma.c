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

#include "fmprb.h"
#include "gamma.h"

void
fmprb_digamma(fmprb_t y, const fmprb_t x, long prec)
{
    int reflect;
    long r, n, wp;
    fmprb_t t, u, v;

    wp = prec + FLINT_BIT_COUNT(prec);

    gamma_stirling_choose_param_fmprb(&reflect, &r, &n, x, 1, 1, wp);

    fmprb_init(t);
    fmprb_init(u);
    fmprb_init(v);

    /* psi(x) = psi((1-x)+r) - h(1-x,r) - pi*cot(pi*x) */
    if (reflect)
    {
        fmprb_sub_ui(t, x, 1, wp);
        fmprb_neg(t, t);
        fmprb_cot_pi(v, x, wp);
        fmprb_const_pi(u, wp);
        fmprb_mul(v, v, u, wp);
        gamma_harmonic_sum_fmprb_ui_bsplit(u, t, r, wp);
        fmprb_add(v, v, u, wp);
        fmprb_add_ui(t, t, r, wp);
        gamma_stirling_eval_fmprb(u, t, n, 1, wp);
        fmprb_sub(y, u, v, wp);
    }
    else
    {
        fmprb_add_ui(t, x, r, wp);
        gamma_stirling_eval_fmprb(u, t, n, 1, wp);
        gamma_harmonic_sum_fmprb_ui_bsplit(t, x, r, wp);
        fmprb_sub(y, u, t, prec);
    }

    fmprb_clear(t);
    fmprb_clear(u);
    fmprb_clear(v);
}

