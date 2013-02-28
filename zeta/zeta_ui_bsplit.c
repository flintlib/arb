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

/* With parameter n, the error is bounded by 3/(3+sqrt(8))^n */
#define ERROR_A 1.5849625007211561815 /* log2(3) */
#define ERROR_B 2.5431066063272239453 /* log2(3+sqrt(8)) */

typedef struct
{
    fmprb_t A;
    fmprb_t B;
    fmprb_t C;
    fmprb_t D;
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
    fmprb_clear(S->Q1);
    fmprb_clear(S->Q2);
    fmprb_clear(S->Q3);
}


static __inline__ void
zeta_coeff_k(zeta_bsplit_t S, long k, long n, long s)
{
    fmprb_set_si(S->D, 2 * (n + k));
    fmprb_mul_si(S->D, S->D, n - k, FMPR_PREC_EXACT);
    fmprb_set_si(S->Q1, k + 1);
    fmprb_mul_si(S->Q1, S->Q1, 2*k + 1, FMPR_PREC_EXACT);

    if (k == 0)
    {
        fmprb_zero(S->A);
        fmprb_one(S->Q2);
    }
    else
    {
        fmprb_set_si(S->A, k % 2 ? 1 : -1);
        fmprb_mul(S->A, S->A, S->Q1, FMPR_PREC_EXACT);
        fmprb_ui_pow_ui(S->Q2, k, s, FMPR_PREC_EXACT);
    }

    fmprb_mul(S->Q3, S->Q1, S->Q2, FMPR_PREC_EXACT);
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

        fmprb_mul(L->B, L->B, R->D, bits);
        fmprb_addmul(L->B, L->A, R->C, bits);

        fmprb_mul(L->B, L->B, R->Q2, bits);
        fmprb_addmul(L->B, R->B, L->Q3, bits);

        fmprb_mul(L->A, L->A, R->Q3, bits);
        fmprb_addmul(L->A, R->A, L->Q3, bits);

        fmprb_mul(L->C, L->C, R->D, bits);
        fmprb_addmul(L->C, R->C, L->Q1, bits);

        if (cont)
        {
            fmprb_mul(L->D, L->D, R->D, bits);
            fmprb_mul(L->Q2, L->Q2, R->Q2, bits);
        }

        fmprb_mul(L->Q1, L->Q1, R->Q1, bits);
        fmprb_mul(L->Q3, L->Q3, R->Q3, bits);

        zeta_bsplit_clear(R);
    }
}

/* The error for eta(s) is bounded by 3/(3+sqrt(8))^n */
void
borwein_error(fmpr_t err, long n)
{
    fmpr_sqrt_ui(err, 8, FMPRB_RAD_PREC, FMPR_RND_DOWN);
    fmpr_add_ui(err, err, 3, FMPRB_RAD_PREC, FMPR_RND_DOWN);
    fmpr_pow_sloppy_ui(err, err, n, FMPRB_RAD_PREC, FMPR_RND_DOWN);
    fmpr_ui_div(err, 3, err, FMPRB_RAD_PREC, FMPR_RND_UP);
}

void
fmprb_zeta_ui_bsplit(fmprb_t x, ulong s, long prec)
{
    zeta_bsplit_t sum;
    fmpr_t err;
    long wp, n;

    /* zeta(0) = -1/2 */
    if (s == 0)
    {
        fmpr_set_si_2exp_si(fmprb_midref(x), -1, -1);
        fmpr_zero(fmprb_radref(x));
        return;
    }

    if (s == 1)
    {
        printf("fmprb_zeta_ui_bsplit: zeta(1)");
        abort();
    }

    n = prec / ERROR_B + 2;
    wp = prec + 30;

    zeta_bsplit_init(sum);
    zeta_bsplit(sum, 0, n + 1, n, s, 0, wp);

    /*  A/Q3 - B/Q3 / (C/Q1) = (A*C - B*Q1) / (Q3*C)    */
    fmprb_mul(sum->A, sum->A, sum->C, wp);
    fmprb_mul(sum->B, sum->B, sum->Q1, wp);
    fmprb_sub(sum->A, sum->A, sum->B, wp);
    fmprb_mul(sum->Q3, sum->Q3, sum->C, wp);
    fmprb_div(sum->C, sum->A, sum->Q3, wp);

    fmpr_init(err);
    borwein_error(err, n);
    fmprb_add_error_fmpr(sum->C, err);
    fmpr_clear(err);

    /* convert from eta(s) to zeta(s) */
    fmprb_div_2expm1_ui(x, sum->C, s - 1, wp);
    fmprb_mul_2exp_si(x, x, s - 1);

    zeta_bsplit_clear(sum);
}

