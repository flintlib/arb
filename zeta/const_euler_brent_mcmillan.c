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

#include "zeta.h"

typedef struct
{
    fmprb_t P;
    fmprb_t Q;
    fmprb_t T;
    fmprb_t C;
    fmprb_t D;
    fmprb_t V;
} euler_bsplit_struct;

typedef euler_bsplit_struct euler_bsplit_t[1];

static void euler_bsplit_init(euler_bsplit_t s)
{
    fmprb_init(s->P);
    fmprb_init(s->Q);
    fmprb_init(s->T);
    fmprb_init(s->C);
    fmprb_init(s->D);
    fmprb_init(s->V);
}

static void euler_bsplit_clear(euler_bsplit_t s)
{
    fmprb_clear(s->P);
    fmprb_clear(s->Q);
    fmprb_clear(s->T);
    fmprb_clear(s->C);
    fmprb_clear(s->D);
    fmprb_clear(s->V);
}

static void
euler_bsplit_1_merge(euler_bsplit_t s, euler_bsplit_t L, euler_bsplit_t R,
                    long wp, int cont)
{
    fmprb_t t, u, v;

    fmprb_init(t);
    fmprb_init(u);
    fmprb_init(v);

    if (cont)
        fmprb_mul(s->P, L->P, R->P, wp);

    fmprb_mul(s->Q, L->Q, R->Q, wp);
    fmprb_mul(s->D, L->D, R->D, wp);

    /* T = LP RT + RQ LT*/
    fmprb_mul(t, L->P, R->T, wp);
    fmprb_mul(v, R->Q, L->T, wp);
    fmprb_add(s->T, t, v, wp);

    /* C = LC RD + RC LD */
    if (cont)
    {
        fmprb_mul(s->C, L->C, R->D, wp);
        fmprb_addmul(s->C, R->C, L->D, wp);
    }

    /* V = RD (RQ LV + LC LP RT) + LD LP RV */
    fmprb_mul(u, L->P, R->V, wp);
    fmprb_mul(u, u, L->D, wp);
    fmprb_mul(v, R->Q, L->V, wp);
    fmprb_addmul(v, t, L->C, wp);
    fmprb_mul(v, v, R->D, wp);
    fmprb_add(s->V, u, v, wp);

    fmprb_clear(t);
    fmprb_clear(u);
    fmprb_clear(v);
}

void
euler_bsplit_1(euler_bsplit_t s, long n1, long n2, long N, long wp, int cont)
{
    if (n2 - n1 == 1)
    {
        fmprb_set_si(s->P, N); /* p = N^2 todo: shift optimization */
        fmprb_mul(s->P, s->P, s->P, wp);
        fmprb_set_si(s->Q, n1 + 1); /* q = (k + 1)^2 */
        fmprb_mul(s->Q, s->Q, s->Q, wp);
        fmprb_set_si(s->C, 1);
        fmprb_set_si(s->D, n1 + 1);
        fmprb_set(s->T, s->P);
        fmprb_set(s->V, s->P);
    }
    else
    {
        euler_bsplit_t L, R;
        long m = (n1 + n2) / 2;

        euler_bsplit_init(L);
        euler_bsplit_init(R);
        euler_bsplit_1(L, n1, m, N, wp, 1);
        euler_bsplit_1(R, m, n2, N, wp, 1);
        euler_bsplit_1_merge(s, L, R, wp, cont);
        euler_bsplit_clear(L);
        euler_bsplit_clear(R);
    }
}

void
euler_bsplit_2(fmprb_t P, fmprb_t Q, fmprb_t T, long n1, long n2,
                        long N, long wp, int cont)
{
    if (n2 - n1 == 1)
    {
        if (n1 == 0)
        {
            fmprb_set_si(P, 1);
            fmprb_set_si(Q, 4 * N);
            fmprb_set_si(T, 1);
        }
        else
        {
            fmprb_si_pow_ui(P, 1 - 2*n1, 3, wp);
            fmprb_neg(P, P);

            fmprb_set_si(Q, 32 * n1);
            fmprb_mul_ui(Q, Q, N, wp);
            fmprb_mul_ui(Q, Q, N, wp);
        }

        fmprb_set(T, P);
    }
    else
    {
        fmprb_t P2, Q2, T2;
        long m = (n1 + n2) / 2;

        fmprb_init(P2);
        fmprb_init(Q2);
        fmprb_init(T2);

        euler_bsplit_2(P, Q, T, n1, m, N, wp, 1);
        euler_bsplit_2(P2, Q2, T2, m, n2, N, wp, 1);

        fmprb_mul(T, T, Q2, wp);
        fmprb_mul(T2, T2, P, wp);
        fmprb_add(T, T, T2, wp);

        if (cont)
            fmprb_mul(P, P, P2, wp);

        fmprb_mul(Q, Q, Q2, wp);

        fmprb_clear(P2);
        fmprb_clear(Q2);
        fmprb_clear(T2);
    }
}

void
fmprb_const_euler_brent_mcmillan(fmprb_t res, long prec)
{
    euler_bsplit_t sum;
    fmprb_t t, u, v, P2, T2, Q2;
    long bits, wp, n, nterms1, nterms2;

    bits = prec + 20;
    n = 0.08665 * bits + 1;

    nterms1 = 4.9706258 * n + 1;
    nterms2 = 2 * n + 1;
    wp = bits + FLINT_BIT_COUNT(n);

    euler_bsplit_init(sum);
    fmprb_init(P2);
    fmprb_init(T2);
    fmprb_init(Q2);
    fmprb_init(t);
    fmprb_init(u);
    fmprb_init(v);

    /* Compute S0 = V / (Q * D), I0 = 1 + T / Q */
    euler_bsplit_1(sum, 0, nterms1, n, wp, 0);

    /* Compute K0 = T2 / Q2 */
    euler_bsplit_2(P2, Q2, T2, 0, nterms2, n, wp, 0);

    /* Compute (S0/I0 + K0/I0^2) = (Q2*(Q+T)*V - D*Q^2*T2)/(D*Q2*(Q+T)^2) */
    fmprb_add(v, sum->Q, sum->T, wp);
    fmprb_mul(t, v, Q2, wp);
    fmprb_mul(u, sum->Q, sum->Q, wp);
    fmprb_mul(u, u, T2, wp);
    fmprb_mul(u, u, sum->D, wp);
    fmprb_mul(sum->V, t, sum->V, wp);
    fmprb_sub(sum->V, sum->V, u, wp);
    fmprb_mul(u, sum->D, t, wp);
    fmprb_mul(u, u, v, wp);
    fmprb_div(t, sum->V, u, wp);

    /* subtract log(n) */
    fmprb_log_ui(u, n, wp);
    fmprb_sub(res, t, u, wp);

    /* TODO: add error term */

    fmprb_clear(P2);
    fmprb_clear(T2);
    fmprb_clear(Q2);
    fmprb_clear(t);
    fmprb_clear(u);
    fmprb_clear(v);

    euler_bsplit_clear(sum);
}
