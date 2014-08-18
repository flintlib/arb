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

#include "arb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("log....");
    fflush(stdout);

    flint_randinit(state);

    /* compare with mpfr */
    for (iter = 0; iter < 100000; iter++)
    {
        arb_t a, b;
        fmpq_t q;
        mpfr_t t;
        long prec = 2 + n_randint(state, 200);

        arb_init(a);
        arb_init(b);
        fmpq_init(q);
        mpfr_init2(t, prec + 300);

        do {
            arb_randtest(a, state, 1 + n_randint(state, 200), 10);
        } while (arb_contains_nonpositive(a));

        arb_randtest(b, state, 1 + n_randint(state, 200), 10);
        arb_get_rand_fmpq(q, state, a, 1 + n_randint(state, 200));

        fmpq_get_mpfr(t, q, MPFR_RNDN);
        mpfr_log(t, t, MPFR_RNDN);

        arb_log(b, a, prec);

        if (!arb_contains_mpfr(b, t))
        {
            printf("FAIL: containment\n\n");
            printf("a = "); arb_print(a); printf("\n\n");
            printf("b = "); arb_print(b); printf("\n\n");
            abort();
        }

        arb_log(a, a, prec);

        if (!arb_equal(a, b))
        {
            printf("FAIL: aliasing\n\n");
            abort();
        }

        arb_clear(a);
        arb_clear(b);
        fmpq_clear(q);
        mpfr_clear(t);
    }

    /* compare with mpfr (higher precision) */
    for (iter = 0; iter < 1000; iter++)
    {
        arb_t a, b;
        fmpq_t q;
        mpfr_t t;
        long prec = 2 + n_randint(state, 6000);

        arb_init(a);
        arb_init(b);
        fmpq_init(q);
        mpfr_init2(t, prec + 7000);

        do {
            arb_randtest(a, state, 1 + n_randint(state, 6000), 10);
        } while (arb_contains_nonpositive(a));

        arb_randtest(b, state, 1 + n_randint(state, 6000), 10);
        arb_get_rand_fmpq(q, state, a, 1 + n_randint(state, 200));

        fmpq_get_mpfr(t, q, MPFR_RNDN);
        mpfr_log(t, t, MPFR_RNDN);

        arb_log(b, a, prec);

        if (!arb_contains_mpfr(b, t))
        {
            printf("FAIL: containment\n\n");
            printf("a = "); arb_print(a); printf("\n\n");
            printf("b = "); arb_print(b); printf("\n\n");
            abort();
        }

        arb_log(a, a, prec);

        if (!arb_equal(a, b))
        {
            printf("FAIL: aliasing\n\n");
            abort();
        }

        arb_clear(a);
        arb_clear(b);
        fmpq_clear(q);
        mpfr_clear(t);
    }

    /* test large numbers */
    for (iter = 0; iter < 10000; iter++)
    {
        arb_t a, b, ab, lab, la, lb, lalb;
        long prec = 2 + n_randint(state, 6000);

        arb_init(a);
        arb_init(b);
        arb_init(ab);
        arb_init(lab);
        arb_init(la);
        arb_init(lb);
        arb_init(lalb);

        arb_randtest(a, state, 1 + n_randint(state, 400), 400);
        arb_randtest(b, state, 1 + n_randint(state, 400), 400);

        arb_log(la, a, prec);
        arb_log(lb, b, prec);
        arb_mul(ab, a, b, prec);
        arb_log(lab, ab, prec);
        arb_add(lalb, la, lb, prec);

        if (!arb_overlaps(lab, lalb))
        {
            printf("FAIL: containment\n\n");
            printf("a = "); arb_print(a); printf("\n\n");
            printf("b = "); arb_print(b); printf("\n\n");
            printf("la = "); arb_print(la); printf("\n\n");
            printf("lb = "); arb_print(lb); printf("\n\n");
            printf("ab = "); arb_print(ab); printf("\n\n");
            printf("lab = "); arb_print(lab); printf("\n\n");
            printf("lalb = "); arb_print(lalb); printf("\n\n");
            abort();
        }

        arb_log(a, a, prec);
        if (!arb_overlaps(a, la))
        {
            printf("FAIL: aliasing\n\n");
            abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(ab);
        arb_clear(lab);
        arb_clear(la);
        arb_clear(lb);
        arb_clear(lalb);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

