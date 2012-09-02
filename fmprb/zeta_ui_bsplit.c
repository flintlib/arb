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

/* With parameter n, the error is bounded by 3/(3+sqrt(8))^n */
#define ERROR_A 1.5849625007211561815 /* log2(3) */
#define ERROR_B 2.5431066063272239453 /* log2(3+sqrt(8)) */

typedef struct
{
    fmprb_t A;
    fmprb_t B;
    fmprb_t C;
    fmprb_t D;
    fmprb_t E;
    fmprb_t Q1;
    fmprb_t Q2;
    fmprb_t Q3;
}
zeta_bsplit_state;

typedef zeta_bsplit_state zeta_bsplit_t[1];

static __inline__ void
zeta_bsplit_init(zeta_bsplit_t S)
{
    fmprb_init(S->A);
    fmprb_init(S->B);
    fmprb_init(S->C);
    fmprb_init(S->D);
    fmprb_init(S->E);
    fmprb_init(S->Q1);
    fmprb_init(S->Q2);
    fmprb_init(S->Q3);
}

static __inline__ void
zeta_bsplit_clear(zeta_bsplit_t S)
{
    fmprb_clear(S->A);
    fmprb_clear(S->B);
    fmprb_clear(S->C);
    fmprb_clear(S->D);
    fmprb_clear(S->E);
    fmprb_clear(S->Q1);
    fmprb_clear(S->Q2);
    fmprb_clear(S->Q3);
}


static __inline__ void
zeta_coeff_k(zeta_bsplit_t S, long k, long n, long s)
{
    if (k + 1 < 0)
    {
        fmprb_set_si(S->D, 1);
        fmprb_set_si(S->Q1, 1);
    }
    else if (k + 1 > n)
    {
        fmprb_zero(S->D);
        fmprb_set_si(S->Q1, 1);
    }
    else
    {
        fmprb_set_si(S->D, 2 * (n + (k + 1) - 1));
        fmprb_mul_si(S->D, S->D, n + 1 - (k + 1), FMPR_PREC_EXACT);
        fmprb_set_si(S->Q1, k + 1);
        fmprb_mul_si(S->Q1, S->Q1, 2*(k + 1) - 1, FMPR_PREC_EXACT);
    }

    if (k - 1 < 0)
    {
        fmprb_zero(S->E);
        fmprb_set_si(S->Q2, 1);
    }
    else if (k - 1 >= n)
    {
        fmprb_set_si(S->E, 1);
        fmprb_set_si(S->Q2, 1);
    }
    else
    {
        fmprb_set_si(S->E, ((k - 1) % 2) ? -1 : 1);
        fmprb_set_si(S->Q2, k);
        fmprb_ui_pow_ui(S->Q2, k, s, FMPR_PREC_EXACT);   /* XXX */
    }

    fmprb_mul(S->Q3, S->Q1, S->Q2, FMPR_PREC_EXACT);     /* XXX */
    fmprb_mul(S->A, S->E, S->Q1, FMPR_PREC_EXACT);
    fmprb_zero(S->B);
    fmprb_set(S->C, S->Q1);
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

        zeta_bsplit_init(R);
        zeta_bsplit(R, a, m, n, s, 1, bits);

        fmprb_mul(L->E, L->E, R->Q2, bits);
        fmprb_addmul(L->E, R->E, L->Q2, bits);

        fmprb_mul(L->B, L->B, R->D, bits);
        fmprb_addmul(L->B, L->A, R->C, bits);

        fmprb_mul(L->B, L->B, R->Q2, bits);
        fmprb_addmul(L->B, R->B, L->Q3, bits);

        if (cont)
        {
            fmprb_mul(L->A, L->A, R->Q3, bits);
            fmprb_addmul(L->A, R->A, L->Q3, bits);
        }

        fmprb_mul(L->C, L->C, R->D, bits);
        fmprb_addmul(L->C, R->C, L->Q1, bits);
        fmprb_mul(L->Q2, L->Q2, R->Q2, bits);

        if (cont)
        {
            fmprb_mul(L->D, L->D, R->D, bits);
            fmprb_mul(L->Q1, L->Q1, R->Q1, bits);
            fmprb_mul(L->Q3, L->Q3, R->Q3, bits);
        }

        zeta_bsplit_clear(R);
    }
}

void
fmprb_zeta_ui_bsplit(fmprb_t x, ulong s, long prec)
{
    zeta_bsplit_t sum;
    long wp, n;
    long i;

    /* zeta(0) = -1/2 */
    if (s == 0)
    {
        /* XXX */
        fmpz_set_si(fmpr_manref(fmprb_midref(x)), -1);
        fmpz_set_si(fmpr_expref(fmprb_midref(x)), -1);
        fmpr_zero(fmprb_radref(x));
        return;
    }

    if (s == 1)
    {
        printf("fmprb_zeta_ui_bsplit: zeta(1)");
        abort();
    }

    /* for s > p, zeta(s) - 1 < 2^(-p) */
    if (s != 2 && s > prec)
    {
        fmprb_set_ui(x, 1UL);

        /* XXX: could make this even smaller when s is extremely large */
        fmpz_one(fmpr_manref(fmprb_radref(x)));
        fmpz_set_si(fmpr_expref(fmprb_radref(x)), -prec);

        return;
    }

    n = prec / ERROR_B + 2;
    wp = prec + 30;

    zeta_bsplit_init(sum);
    zeta_bsplit(sum, 0, n + 1, n, s, 0, wp);

    fmprb_mul(sum->E, sum->E, sum->C, wp);
    fmprb_sub(sum->E, sum->E, sum->B, wp);
    fmprb_mul(sum->Q2, sum->Q2, sum->C, wp);
    fmprb_div(sum->C, sum->E, sum->Q2, wp);

    /* D = 1/(1 - 2^(1-s)) */
    {
        fmpz_t t;
        fmpz_init(t);
        for (i = wp; i >= 0; i -= (s - 1))
            fmpz_setbit(t, i);

        fmprb_set_fmpz(sum->D, t);
        fmpz_sub_ui(fmpr_expref(fmprb_midref(sum->D)),
            fmpr_expref(fmprb_midref(sum->D)), wp);

        fmprb_add_error_2exp_si(sum->D, -wp);

        fmpz_clear(t);
    }

    fmprb_mul(x, sum->C, sum->D, wp);

    /* The truncation error is bounded by 3/(3+sqrt(8))^n */
    fmprb_add_error_2exp_si(x, (long) (ERROR_A - ERROR_B * n + 1));

    zeta_bsplit_clear(sum);
}
