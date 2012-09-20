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
#include "fmprb_poly.h"

static void
bsplit(fmprb_t P, fmprb_t Q, fmprb_t B, fmprb_t T, long n1, long n2,
    const fmprb_t anum, const fmprb_t aden, long n, long m, long prec)
{
    if (n2 - n1 == 1)
    {
        long k = n1;

        if (k == 0)
            fmprb_pow_ui(P, aden, m + 1, prec);
        else
            fmprb_set_ui(P, n);

        fmprb_set_ui(Q, (k == 0) ? 1 : k);
        fmprb_set(B, anum);
        fmprb_addmul_ui(B, aden, k, prec);
        fmprb_pow_ui(B, B, m + 1, prec);
        fmprb_mul_si(T, P, (k % 2) ? -1 : 1, prec);
    }
    else
    {
        long mid = (n1 + n2) / 2;

        fmprb_t P2, Q2, B2, T2;

        fmprb_init(P2);
        fmprb_init(Q2);
        fmprb_init(B2);
        fmprb_init(T2);

        bsplit(P,  Q,  B,  T,  n1, mid, anum, aden, n, m, prec);
        bsplit(P2, Q2, B2, T2, mid, n2, anum, aden, n, m, prec);

        if (!fmprb_is_one(B2))
            fmprb_mul(T, T, B2, prec);

        fmprb_mul(T, T, Q2, prec);

        if (!fmprb_is_one(B))
            fmprb_mul(T2, T2, B, prec);

        fmprb_mul(T2, T2, P, prec);
        fmprb_add(T, T, T2, prec);

        fmprb_mul(P, P, P2, prec);
        fmprb_mul(Q, Q, Q2, prec);
        fmprb_mul(B, B, B2, prec);

        fmprb_clear(P2);
        fmprb_clear(Q2);
        fmprb_clear(B2);
        fmprb_clear(T2);
    }
}

static void
karatsuba_error(fmpr_t bound, long s, long n, long r)
{
    fmpr_t t, u, A, B;
    long wp = FMPRB_RAD_PREC;

    if (r < n || s < 2)
        abort();

    fmpr_init(t);
    fmpr_init(u);
    fmpr_init(A);
    fmpr_init(B);

    /* t = log(n)^s */
    fmpr_set_ui(t, n);
    fmpr_log(t, t, wp, FMPR_RND_UP);
    fmpr_pow_sloppy_ui(t, t, s, wp, FMPR_RND_UP);

    /* A = 5/3 * exp(-n) * log(n) ** s */
    fmpr_set_si(u, -n);
    fmpr_exp(u, u, wp, FMPR_RND_UP);
    fmpr_mul(u, u, t, wp, FMPR_RND_UP);
    fmpr_mul_ui(u, u, 5, wp, FMPR_RND_UP);
    fmpr_div_ui(A, u, 3, wp, FMPR_RND_UP);

    /* B = (e/(r+2))^(r+2) * (1 + n^(r+2) * log(n)**2) */
    fmpr_set_ui(u, n);
    fmpr_pow_sloppy_ui(u, u, r + 2, wp, FMPR_RND_UP);
    fmpr_mul(t, t, u, wp, FMPR_RND_UP);
    fmpr_add_ui(t, t, 1UL, wp, FMPR_RND_UP);
    fmpr_set_ui(u, 1UL);
    fmpr_exp(u, u, wp, FMPR_RND_UP);
    fmpr_div_ui(u, u, r + 2, wp, FMPR_RND_UP);
    fmpr_pow_sloppy_ui(u, u, r + 2, wp, FMPR_RND_UP);
    fmpr_mul(B, t, u, wp, FMPR_RND_UP);

    /* A + B */
    fmpr_add(bound, A, B, wp, FMPR_RND_UP);

    fmpr_clear(t);
    fmpr_clear(u);
    fmpr_clear(A);
    fmpr_clear(B);
}

void
fmprb_gamma_fmpq_karatsuba(fmprb_struct * v, const fmpq_t a, long num, long prec)
{
    long s, jmax, m, n, r, wp, bsplit_wp;
    fmprb_struct * Sm, * logs;
    fmprb_t t, u, anum, aden;
    fmprb_t P, Q, B, T;
    fmpr_t truncation_error;

    fmpr_init(truncation_error);

    /* Select parameters (heuristic) */
    jmax = num - 1;
    s = FLINT_MAX(2, jmax);
    wp = prec * 1.01 + 10 + 2*jmax;
    n = FLINT_MAX(wp * 0.693147180559945, 2*s*FLINT_BIT_COUNT(2*s) + 1);
    /* round n to have many trailing zeros, speeding up arithmetic */
    if (n > 256)
    {
        int b = FLINT_BIT_COUNT(n);
        n = ((n >> (b-4)) + 1) << (b-4);
    }
    r = 3.59112147666862*n;

    /* the bound computation is rigorous */
    karatsuba_error(truncation_error, s, n, r);

    /* Heuristic: when s is small, we can truncate in the binary splitting.
       The factor 2.2 is not generally valid and should be replaced
       based on theory. */
    if (15 * s < wp)
        bsplit_wp = 2.2 * wp;
    else
        bsplit_wp = FMPR_PREC_EXACT;

#if 0
    printf("s = %ld n = %ld r = %ld err = ", s, n, r);
        fmpr_printd(truncation_error, 10); printf("\n");
#endif

    fmprb_init(t);
    fmprb_init(u);
    fmprb_init(anum);
    fmprb_init(aden);
    fmprb_init(P);
    fmprb_init(Q);
    fmprb_init(B);
    fmprb_init(T);

    fmprb_set_fmpz(anum, fmpq_numref(a));
    fmprb_set_fmpz(aden, fmpq_denref(a));

    Sm = _fmprb_vec_init(jmax + 1);
    logs = _fmprb_vec_init(jmax + 1);

    /* Compute S sums */
    for (m = 0; m <= jmax; m++) {
        /* 2.2 * wp */
        bsplit(P, Q, B, T, 0, r + 1, anum, aden, n, m, bsplit_wp);
        fmprb_mul(Q, Q, B, wp);
        fmprb_div(Sm + m, T, Q, wp);
#if 0
        printf("%ld %ld ", m, jmax); fmprb_printd(Sm + m, 10); printf("\n");
#endif
        if (m % 2)
            fmprb_neg(Sm + m, Sm + m);
    }

    /* Compute (log(n))^m / m! */
    for (m = 0; m <= jmax; m++)
    {
        if (m == 0) {
            fmprb_one(logs + m);
        } else if (m == 1) {
            fmprb_log_ui(logs + m, n, wp);
        } else {
            fmprb_mul(logs + m, logs + m - 1, logs + 1, wp);
            fmprb_div_ui(logs + m, logs + m, m, wp);
        }
    }

    /* Convolution */
    _fmprb_poly_mullow(v, logs, jmax + 1, Sm, jmax + 1, jmax + 1, wp);

    /* multiply by n^a */
    fmprb_set_ui(t, n);
    fmprb_pow_fmpq(t, t, a, wp);
    for (m = 0; m <= jmax; m++)
        fmprb_mul(v + m, v + m, t, wp);

    /* Add error. TODO: the bound is for the derivatives directly,
       so we should probably divide it by m! ? */
    for (m = 0; m <= jmax; m++)
        fmprb_add_error_fmpr(v + m, truncation_error);

    fmprb_clear(t);
    fmprb_clear(u);
    fmprb_clear(anum);
    fmprb_clear(aden);
    fmprb_clear(P);
    fmprb_clear(Q);
    fmprb_clear(B);
    fmprb_clear(T);
    fmpr_clear(truncation_error);

    _fmprb_vec_clear(Sm, jmax + 1);
    _fmprb_vec_clear(logs, jmax + 1);
}
