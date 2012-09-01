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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "fmprb.h"

#define CONST_A 13591409UL
#define CONST_B 545140134UL
#define CONST_C 640320UL
#define CONST_D 12UL
#define BITS_PER_TERM 47.1104131382158 /* log2(640320^3 / (2^6 * 3^3)) */

void
chudnovsky_bsplit(fmprb_t G, fmprb_t P, fmprb_t Q, long a, long b, long wp, int cont)
{
    if (b - a == 1)
    {
        /* g = (6*b-5)*(2*b-1)*(6*b-1) */
        fmprb_set_si(G, 6*b - 5);
        fmprb_mul_si(G, G, 2*b-1, FMPR_PREC_EXACT);
        fmprb_mul_si(G, G, 6*b-1, FMPR_PREC_EXACT);

        /* p = C^3 * b^3 / 24 */
        fmprb_set_si(P, CONST_C * CONST_C * CONST_C / 24); /* xxx: fix for 32-bit */
        fmprb_mul_ui(P, P, b, FMPR_PREC_EXACT);
        fmprb_mul_ui(P, P, b, FMPR_PREC_EXACT);
        fmprb_mul_ui(P, P, b, FMPR_PREC_EXACT);

        /* (-1)^b * g * (A + B*b) */
        {
            fmpz_t t;
            fmpz_init(t);
            fmpz_set_ui(t, CONST_B);
            fmpz_mul_ui(t, t, b);
            fmpz_add_ui(t, t, CONST_A);
            fmprb_mul_fmpz(Q, G, t, FMPR_PREC_EXACT);
            fmpz_clear(t);
        }

        if (b % 2)
            fmprb_neg(Q, Q);
    }
    else
    {
        long m;
        fmprb_t G2, P2, Q2;

        m = (a + b) / 2;

        fmprb_init(G2);
        fmprb_init(P2);
        fmprb_init(Q2);

        chudnovsky_bsplit(G, P, Q, a, m, wp, 1);
        chudnovsky_bsplit(G2, P2, Q2, m, b, wp, 1);

        fmprb_mul(Q, Q, P2, wp);
        fmprb_addmul(Q, Q2, G, wp);

        fmprb_mul(P, P, P2, wp);

        if (cont)
            fmprb_mul(G, G, G2, wp);

        fmprb_clear(G2);
        fmprb_clear(P2);
        fmprb_clear(Q2);
    }
}

void
fmprb_const_pi_chudnovsky(fmprb_t x, long prec)
{
    long wp, N;
    fmprb_t G, P, Q;

    wp = prec + 32;

    N = wp / BITS_PER_TERM + 1;

    fmprb_init(G);
    fmprb_init(P);
    fmprb_init(Q);

    chudnovsky_bsplit(G, P, Q, 0, N, wp, 0);

    fmprb_sqrt_ui(G, CONST_C, wp);
    fmprb_mul_si(G, G, CONST_C, wp);
    fmprb_mul(G, G, P, wp);
    fmprb_mul_si(P, P, CONST_A, wp);
    fmprb_add(Q, Q, P, wp);
    fmprb_mul_si(Q, Q, CONST_D, wp);

    fmprb_div(x, G, Q, wp);

    /* TODO: add error term */

    fmprb_clear(G);
    fmprb_clear(P);
    fmprb_clear(Q);
}