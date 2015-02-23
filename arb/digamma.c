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

#include "arb.h"
#include "arith.h"

void arb_gamma_stirling_choose_param(int * reflect, long * r, long * n,
    const arb_t x, int use_reflect, int digamma, long prec);

void arb_gamma_stirling_eval(arb_t s, const arb_t z, long nterms, int digamma, long prec);

void
arb_digamma(arb_t y, const arb_t x, long prec)
{
    int reflect;
    long r, n, wp;
    arb_t t, u, v;

    if (arb_is_exact(x))
    {
        const arf_struct * mid = arb_midref(x);

        if (arf_is_special(mid))
        {
            arb_indeterminate(y);
            return;
        }
        else if (arf_is_int(mid))
        {
            if (arf_sgn(mid) < 0)
            {
                arb_indeterminate(y);
                return;
            }
            else if (arf_cmpabs_ui(mid, 30 + prec / 2) < 0)
            {
                fmpq_t h;

                arb_init(t);
                fmpq_init(h);

                /* should change to fmpq_harmonic_ui when flint upgrades */
                arith_harmonic_number(h, arf_get_si(mid, ARF_RND_DOWN) - 1);
                arb_set_fmpq(y, h, prec + 2);
                arb_const_euler(t, prec + 2);
                arb_sub(y, y, t, prec);

                arb_clear(t);
                fmpq_clear(h);

                return;
            }
        }
    }

    wp = prec + FLINT_BIT_COUNT(prec);

    arb_gamma_stirling_choose_param(&reflect, &r, &n, x, 1, 1, wp);

    arb_init(t);
    arb_init(u);
    arb_init(v);

    /* psi(x) = psi((1-x)+r) - h(1-x,r) - pi*cot(pi*x) */
    if (reflect)
    {
        arb_sub_ui(t, x, 1, wp);
        arb_neg(t, t);
        arb_cot_pi(v, x, wp);
        arb_const_pi(u, wp);
        arb_mul(v, v, u, wp);
        arb_rising2_ui(y, u, t, r, wp);
        arb_div(u, u, y, wp);
        arb_add(v, v, u, wp);
        arb_add_ui(t, t, r, wp);
        arb_gamma_stirling_eval(u, t, n, 1, wp);
        arb_sub(y, u, v, wp);
    }
    else
    {
        arb_add_ui(t, x, r, wp);
        arb_gamma_stirling_eval(u, t, n, 1, wp);
        arb_rising2_ui(y, t, x, r, wp);
        arb_div(t, t, y, wp);
        arb_sub(y, u, t, prec);
    }

    arb_clear(t);
    arb_clear(u);
    arb_clear(v);
}

