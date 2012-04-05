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

#define CONST_A 13591409UL
#define CONST_B 545140134UL
#define CONST_C 640320UL
#define CONST_D 12UL
#define BITS_PER_TERM 47.1104131382158  /* log2(640320^3 / (2^6 * 3^3)) */

void
chudnovsky_bsplit(arb_t G, arb_t P, arb_t Q, long a, long b, long wp, int cont)
{
    if (b - a == 1)
    {
        /* g = (6*b-5)*(2*b-1)*(6*b-1) */
        arb_set_si(G, 6*b - 5);
        fmpz_mul_si(arb_midref(G), arb_midref(G), 2*b-1);
        fmpz_mul_si(arb_midref(G), arb_midref(G), 6*b-1);

        /* p = C^3 * b^3 / 24 */
        arb_set_si(P, CONST_C * CONST_C * CONST_C / 24);  /* xxx: fix for 32-bit */
        fmpz_mul_ui(arb_midref(P), arb_midref(P), b);
        fmpz_mul2_uiui(arb_midref(P), arb_midref(P), b, b);

        /* (-1)^b * g * (A + B*b) */
        arb_set_si(Q, CONST_A);
        fmpz_addmul_ui(arb_midref(Q), (fmpz *) &b, CONST_B);
        arb_mul(Q, Q, G);
        if (b % 2)
            fmpz_neg(arb_midref(Q), arb_midref(Q));
    }
    else
    {
        long m;
        arb_t G2, P2, Q2;

        m = (a + b) / 2;

        arb_init(G2, wp);
        arb_init(P2, wp);
        arb_init(Q2, wp);

        chudnovsky_bsplit(G, P, Q, a, m, wp, 1);
        chudnovsky_bsplit(G2, P2, Q2, m, b, wp, 1);

        arb_mul(Q, Q, P2);
        arb_addmul(Q, Q2, G);

        arb_mul(P, P, P2);

        if (cont)
            arb_mul(G, G, G2);

        arb_clear(G2);
        arb_clear(P2);
        arb_clear(Q2);
    }
}

void
arb_const_pi_chudnovsky(arb_t x)
{
    long prec, wp, N;
    arb_t G, P, Q;

    prec = arb_prec(x);
    wp = prec + 32;

    N = wp / BITS_PER_TERM + 1;

    arb_init(G, wp);
    arb_init(P, wp);
    arb_init(Q, wp);

    chudnovsky_bsplit(G, P, Q, 0, N, wp, 0);

    arb_sqrt_ui(G, CONST_C);
    arb_mul_si(G, G, CONST_C);
    arb_mul(G, G, P);
    arb_mul_si(P, P, CONST_A);
    arb_add(Q, Q, P);
    arb_mul_si(Q, Q, CONST_D);
    arb_div(x, G, Q);

    /* TODO: add error term */

    arb_clear(G);
    arb_clear(P);
    arb_clear(Q);
}
