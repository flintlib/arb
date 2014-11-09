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

#include "acb.h"

void acb_gamma_stirling_choose_param(int * reflect, long * r, long * n,
    const acb_t x, int use_reflect, int digamma, long prec);

void acb_gamma_stirling_eval(acb_t s, const acb_t z, long nterms, int digamma, long prec);

void
acb_digamma(acb_t y, const acb_t x, long prec)
{
    int reflect;
    long r, n, wp;
    acb_t t, u, v;

    if (acb_is_real(x))
    {
        arb_digamma(acb_realref(y), acb_realref(x), prec);
        arb_zero(acb_imagref(y));
        return;
    }

    wp = prec + FLINT_BIT_COUNT(prec);

    acb_gamma_stirling_choose_param(&reflect, &r, &n, x, 1, 1, wp);

    acb_init(t);
    acb_init(u);
    acb_init(v);

    /* psi(x) = psi((1-x)+r) - h(1-x,r) - pi*cot(pi*x) */
    if (reflect)
    {
        acb_sub_ui(t, x, 1, wp);
        acb_neg(t, t);
        acb_cot_pi(v, x, wp);
        arb_const_pi(acb_realref(u), wp);
        acb_mul_arb(v, v, acb_realref(u), wp);
        acb_rising2_ui(y, u, t, r, wp);
        acb_div(u, u, y, wp);
        acb_add(v, v, u, wp);
        acb_add_ui(t, t, r, wp);
        acb_gamma_stirling_eval(u, t, n, 1, wp);
        acb_sub(y, u, v, wp);
    }
    else
    {
        acb_add_ui(t, x, r, wp);
        acb_gamma_stirling_eval(u, t, n, 1, wp);
        acb_rising2_ui(y, t, x, r, wp);
        acb_div(t, t, y, wp);
        acb_sub(y, u, t, prec);
    }

    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
}

