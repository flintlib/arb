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

typedef struct
{
    arb_t P;
    arb_t Q;
    arb_t T;
    arb_t C;
    arb_t D;
    arb_t V;
} euler_bsplit_struct;

typedef euler_bsplit_struct euler_bsplit_t[1];

static void euler_bsplit_init(euler_bsplit_t s, long wp)
{
    arb_init(s->P, wp);
    arb_init(s->Q, wp);
    arb_init(s->T, wp);
    arb_init(s->C, wp);
    arb_init(s->D, wp);
    arb_init(s->V, wp);
}

static void euler_bsplit_clear(euler_bsplit_t s)
{
    arb_clear(s->P);
    arb_clear(s->Q);
    arb_clear(s->T);
    arb_clear(s->C);
    arb_clear(s->D);
    arb_clear(s->V);
}

static void
euler_bsplit_1_merge(euler_bsplit_t s, euler_bsplit_t L, euler_bsplit_t R,
                    long wp, int cont)
{
    arb_t t, u, v;

    arb_init(t, wp);
    arb_init(u, wp);
    arb_init(v, wp);

    if (cont)
        arb_mul(s->P, L->P, R->P);

    arb_mul(s->Q, L->Q, R->Q);
    arb_mul(s->D, L->D, R->D);

    /* T = LP RT + RQ LT*/
    arb_mul(t, L->P, R->T);
    arb_mul(v, R->Q, L->T);
    arb_add(s->T, t, v);

    /* C = LC RD + RC LD */
    if (cont)
    {
        arb_mul(s->C, L->C, R->D);
        arb_addmul(s->C, R->C, L->D);
    }

    /* V = RD (RQ LV + LC LP RT) + LD LP RV */
    arb_mul(u, L->P, R->V);
    arb_mul(u, u, L->D);
    arb_mul(v, R->Q, L->V);
    arb_addmul(v, t, L->C);
    arb_mul(v, v, R->D);
    arb_add(s->V, u, v);

    arb_clear(t);
    arb_clear(u);
    arb_clear(v);
}

void
euler_bsplit_1(euler_bsplit_t s, long n1, long n2, long N, long wp, int cont)
{
    if (n2 - n1 == 1)
    {
        arb_set_si(s->P, N);        /* p = N^2  todo: shift optimization */
        arb_mul(s->P, s->P, s->P);
        arb_set_si(s->Q, n1 + 1);   /* q = (k + 1)^2 */
        arb_mul(s->Q, s->Q, s->Q);
        arb_set_si(s->C, 1);
        arb_set_si(s->D, n1 + 1);
        arb_set(s->T, s->P);
        arb_set(s->V, s->P);
    }
    else
    {
        euler_bsplit_t L, R;
        long m = (n1 + n2) / 2;

        euler_bsplit_init(L, wp);
        euler_bsplit_init(R, wp);
        euler_bsplit_1(L, n1, m, N, wp, 1);
        euler_bsplit_1(R, m, n2, N, wp, 1);
        euler_bsplit_1_merge(s, L, R, wp, cont);
        euler_bsplit_clear(L);
        euler_bsplit_clear(R);
    }
}

void
euler_bsplit_2(arb_t P, arb_t Q, arb_t T, long n1, long n2,
                        long N, long wp, int cont)
{
    if (n2 - n1 == 1)
    {
        if (n1 == 0)
        {
            arb_set_si(P, 1);
            arb_set_si(Q, 4 * N);
            arb_set_si(T, 1);
        }
        else
        {
            arb_zero(P);
            fmpz_set_si(arb_midref(P), 1 - 2*n1);
            fmpz_pow_ui(arb_midref(P), arb_midref(P), 3);
            fmpz_neg(arb_midref(P), arb_midref(P));
            arb_zero(Q);
            fmpz_set_si(arb_midref(Q), 32 * n1);
            fmpz_mul_ui(arb_midref(Q), arb_midref(Q), N);
            fmpz_mul_ui(arb_midref(Q), arb_midref(Q), N);
        }

        arb_set(T, P);
    }
    else
    {
        arb_t P2, Q2, T2;
        long m = (n1 + n2) / 2;

        arb_init(P2, wp);
        arb_init(Q2, wp);
        arb_init(T2, wp);

        euler_bsplit_2(P,  Q,  T,  n1, m, N, wp, 1);
        euler_bsplit_2(P2, Q2, T2, m, n2, N, wp, 1);

        arb_mul(T, T, Q2);
        arb_mul(T2, T2, P);
        arb_add(T, T, T2);

        if (cont)
            arb_mul(P, P, P2);

        arb_mul(Q, Q, Q2);

        arb_clear(P2);
        arb_clear(Q2);
        arb_clear(T2);
    }
}

void
arb_const_euler_brent_mcmillan(arb_t res)
{
    euler_bsplit_t sum;
    arb_t t, u, v, P2, T2, Q2;
    long bits, wp, n, nterms1, nterms2;

    bits = arb_prec(res) + 20;
    n = 0.08665 * bits + 1;

    nterms1 = 4.9706258 * n + 1;
    nterms2 = 2 * n + 1;
    wp = bits + FLINT_BIT_COUNT(n);

    euler_bsplit_init(sum, wp);
    arb_init(P2, wp);
    arb_init(T2, wp);
    arb_init(Q2, wp);
    arb_init(t, wp);
    arb_init(u, wp);
    arb_init(v, wp);

    /* Compute S0 = V / (Q * D), I0 = 1 + T / Q */
    euler_bsplit_1(sum, 0, nterms1, n, wp, 0);

    /* Compute K0 = T2 / Q2 */
    euler_bsplit_2(P2, Q2, T2, 0, nterms2, n, wp, 0);

    /* Compute (S0/I0 + K0/I0^2) = (Q2*(Q+T)*V - D*Q^2*T2)/(D*Q2*(Q+T)^2) */
    arb_add(v, sum->Q, sum->T);
    arb_mul(t, v, Q2);
    arb_mul(u, sum->Q, sum->Q);
    arb_mul(u, u, T2);
    arb_mul(u, u, sum->D);
    arb_mul(sum->V, t, sum->V);
    arb_sub(sum->V, sum->V, u);
    arb_mul(u, sum->D, t);
    arb_mul(u, u, v);
    arb_div(t, sum->V, u);

    /* subtract log(n) */
    arb_log_ui(u, n);
    arb_sub(res, t, u);

    /* TODO: add error term */

    arb_clear(P2);
    arb_clear(T2);
    arb_clear(Q2);
    arb_clear(t);
    arb_clear(u);
    arb_clear(v);

    euler_bsplit_clear(sum);
}
