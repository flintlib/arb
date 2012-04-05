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

/* With parameter n, the error is bounded by 3/(3+sqrt(8))^n */
#define ERROR_A 1.5849625007211561815  /* log2(3) */
#define ERROR_B 2.5431066063272239453  /* log2(3+sqrt(8)) */

typedef struct
{
    arb_t A;
    arb_t B;
    arb_t C;
    arb_t D;
    arb_t E;
    arb_t Q1;
    arb_t Q2;
    arb_t Q3;
}
zeta_bsplit_state;

typedef zeta_bsplit_state zeta_bsplit_t[1];

static __inline__ void
zeta_bsplit_init(zeta_bsplit_t S, long prec)
{
    arb_init(S->A, prec);
    arb_init(S->B, prec);
    arb_init(S->C, prec);
    arb_init(S->D, prec);
    arb_init(S->E, prec);
    arb_init(S->Q1, prec);
    arb_init(S->Q2, prec);
    arb_init(S->Q3, prec);
}

static __inline__ void
zeta_bsplit_clear(zeta_bsplit_t S)
{
    arb_clear(S->A);
    arb_clear(S->B);
    arb_clear(S->C);
    arb_clear(S->D);
    arb_clear(S->E);
    arb_clear(S->Q1);
    arb_clear(S->Q2);
    arb_clear(S->Q3);
}


static __inline__ void
zeta_coeff_k(zeta_bsplit_t S, long k, long n, long s)
{
    if (k + 1 < 0)
    {
        arb_set_si(S->D, 1);
        arb_set_si(S->Q1, 1);
    }
    else if (k + 1 > n)
    {
        arb_zero(S->D);
        arb_set_si(S->Q1, 1);
    }
    else
    {
        arb_set_si(S->D, 2 * (n + (k + 1) - 1));
        arb_mul_si(S->D, S->D, n + 1 - (k + 1));
        arb_set_si(S->Q1, k + 1);
        arb_mul_si(S->Q1, S->Q1, 2*(k + 1) - 1);
    }

    if (k - 1 < 0)
    {
        arb_zero(S->E);
        arb_set_si(S->Q2, 1);
    }
    else if (k - 1 >= n)
    {
        arb_set_si(S->E, 1);
        arb_set_si(S->Q2, 1);
    }
    else
    {
        arb_set_si(S->E, ((k - 1) % 2) ? -1 : 1);
        arb_set_si(S->Q2, k);
        fmpz_pow_ui(arb_midref(S->Q2), arb_midref(S->Q2), s); /* XXX */
    }

    arb_mul(S->Q3, S->Q1, S->Q2);
    arb_mul(S->A, S->E, S->Q1);
    arb_zero(S->B);
    arb_set(S->C, S->Q1);
}

static void
zeta_bsplit(zeta_bsplit_t L, long a, long b,
    long n, long s, int cont, long bits)
{
    if (a + 1 == b)
    {
        zeta_coeff_k(L, a, n, s);
    }
    else
    {
        zeta_bsplit_t R;

        long m = (a + b) / 2;

        zeta_bsplit(L, m, b, n, s, 1, bits);

        zeta_bsplit_init(R, bits);
        zeta_bsplit(R, a, m, n, s, 1, bits);

        arb_mul(L->E, L->E, R->Q2);
        arb_addmul(L->E, R->E, L->Q2);

        arb_mul(L->B, L->B, R->D);
        arb_addmul(L->B, L->A, R->C);

        arb_mul(L->B, L->B, R->Q2);
        arb_addmul(L->B, R->B, L->Q3);

        if (cont)
        {
            arb_mul(L->A, L->A, R->Q3);
            arb_addmul(L->A, R->A, L->Q3);
        }

        arb_mul(L->C, L->C, R->D);
        arb_addmul(L->C, R->C, L->Q1);
        arb_mul(L->Q2, L->Q2, R->Q2);

        if (cont)
        {
            arb_mul(L->D, L->D, R->D);
            arb_mul(L->Q1, L->Q1, R->Q1);
            arb_mul(L->Q3, L->Q3, R->Q3);
        }

        zeta_bsplit_clear(R);
    }
}

void
arb_zeta_ui_bsplit(arb_t x, ulong s)
{
    zeta_bsplit_t sum;
    long prec, wp, n;
    long i;

    /* zeta(0) = -1/2 */
    if (s == 0)
    {
        fmpz_set_si(arb_midref(x), -1);
        fmpz_set_si(arb_expref(x), -1);
        fmpz_zero(arb_radref(x));
        return;
    }

    if (s == 1)
    {
        printf("arb_zeta_ui_bsplit: zeta(1)");
        abort();
    }

    prec = arb_prec(x);

    /* for s > p, zeta(s) - 1 < 2^(-p) */
    if (s != 2 && s > prec)
    {
        fmpz_one(arb_midref(x));
        fmpz_mul_2exp(arb_midref(x), arb_midref(x), prec);
        fmpz_set_si(arb_expref(x), -prec);
        fmpz_one(arb_radref(x));
        return;
    }

    n = prec / ERROR_B + 2;
    wp = prec + 30;

    zeta_bsplit_init(sum, wp);
    zeta_bsplit(sum, 0, n + 1, n, s, 0, wp);

    arb_mul(sum->E, sum->E, sum->C);
    arb_sub(sum->E, sum->E, sum->B);
    arb_mul(sum->Q2, sum->Q2, sum->C);
    arb_div(sum->C, sum->E, sum->Q2);

    /* D = 1/(1 - 2^(1-s)) */
    fmpz_zero(arb_midref(sum->D));
    for (i = wp; i >= 0; i -= (s - 1))
        fmpz_setbit(arb_midref(sum->D), i);
    fmpz_set_si(arb_expref(sum->D), -wp);
    fmpz_set_ui(arb_radref(sum->D), 1UL);

    arb_mul(x, sum->C, sum->D);

    /* The error is bounded by 3/(3+sqrt(8))^n */
    arb_add_error_2exp(x, (long) (ERROR_A - ERROR_B * n + 1));

    zeta_bsplit_clear(sum);
}
