/*=============================================================================

    This file is part of ARB.

    ARB is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    ARB is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ARB; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "arb.h"

static void
zeta3_bsplit(arb_t P, arb_t Q, arb_t T, long a, long b, long wp, int cont)
{
    if (b - a == 1)
    {
        arb_set_si(P, a ? a : 1);
        fmpz_pow_ui(arb_midref(P), arb_midref(P), 5);

        arb_set_si(Q, 2*a + 1);
        fmpz_pow_ui(arb_midref(Q), arb_midref(Q), 5);
        fmpz_mul_ui(arb_midref(Q), arb_midref(Q), a ? 32 : 64);

        arb_set_si(T, 205*a*a + 250*a + 77);     /* xxx: overflow */
        if (a % 2)
            fmpz_neg(arb_midref(T), arb_midref(T));
        arb_mul(T, T, P);
    }
    else
    {
        long m;
        arb_t P2, Q2, T2;

        m = (a + b) / 2;

        arb_init(P2, wp);
        arb_init(Q2, wp);
        arb_init(T2, wp);

        zeta3_bsplit(P, Q, T, a, m, wp, 1);
        zeta3_bsplit(P2, Q2, T2, m, b, wp, 1);

        arb_mul(T, T, Q2);
        arb_addmul(T, P, T2);
        arb_mul(Q, Q, Q2);
        if (cont)
            arb_mul(P, P, P2);

        arb_clear(P2);
        arb_clear(Q2);
        arb_clear(T2);
    }
}

void
arb_const_zeta3_bsplit(arb_t x)
{
    long prec, wp, N;
    arb_t P, Q, T;

    prec = arb_prec(x);
    wp = prec + 10 + FLINT_BIT_COUNT(prec);
    N = prec / 10 + 2;

    arb_init(P, wp);
    arb_init(Q, wp);
    arb_init(T, wp);

    zeta3_bsplit(P, Q, T, 0, N, wp, 0);
    arb_div(x, T, Q);

    arb_add_error_2exp(x, -10 * (N - 1));

    arb_clear(P);
    arb_clear(Q);
    arb_clear(T);
}
