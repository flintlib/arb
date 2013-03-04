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
        fmprb_t a, b;
        fmpq_t q;
        mpfr_t t;
        long prec = 2 + n_randint(state, 200);

        fmprb_init(a);
        fmprb_init(b);
        fmpq_init(q);
        mpfr_init2(t, prec + 200);

        do {
            fmprb_randtest(a, state, 1 + n_randint(state, 200), 10);
        } while (fmprb_contains_nonpositive(a));

        fmprb_randtest(b, state, 1 + n_randint(state, 200), 10);
        fmprb_get_rand_fmpq(q, state, a, 1 + n_randint(state, 200));

        fmpq_get_mpfr(t, q, MPFR_RNDN);
        mpfr_log(t, t, MPFR_RNDN);

        fmprb_log(b, a, prec);

        if (!fmprb_contains_mpfr(b, t))
        {
            printf("FAIL: containment\n\n");
            printf("a = "); fmprb_print(a); printf("\n\n");
            printf("b = "); fmprb_print(b); printf("\n\n");
            abort();
        }

        fmprb_log(a, a, prec);

        if (!fmprb_equal(a, b))
        {
            printf("FAIL: aliasing\n\n");
            abort();
        }

        fmprb_clear(a);
        fmprb_clear(b);
        fmpq_clear(q);
        mpfr_clear(t);
    }

    /* test large numbers */
    for (iter = 0; iter < 10000; iter++)
    {
        fmprb_t a, b, ab, lab, la, lb, lalb;
        long prec = 2 + n_randint(state, 400);

        fmprb_init(a);
        fmprb_init(ab);
        fmprb_init(lab);
        fmprb_init(la);
        fmprb_init(lb);
        fmprb_init(lalb);

        do {
            fmprb_randtest(a, state, 1 + n_randint(state, 400), 400);
        } while (fmprb_contains_nonpositive(a));

        do {
            fmprb_randtest(b, state, 1 + n_randint(state, 400), 400);
        } while (fmprb_contains_nonpositive(a));

        fmprb_log(la, a, prec);
        fmprb_log(lb, b, prec);
        fmprb_mul(ab, a, b, prec);
        fmprb_log(lab, ab, prec);
        fmprb_add(lalb, la, lb, prec);

        if (!fmprb_overlaps(lab, lalb))
        {
            printf("FAIL: containment\n\n");
            printf("a = "); fmprb_print(a); printf("\n\n");
            printf("b = "); fmprb_print(b); printf("\n\n");
            printf("la = "); fmprb_print(la); printf("\n\n");
            printf("lb = "); fmprb_print(lb); printf("\n\n");
            printf("ab = "); fmprb_print(ab); printf("\n\n");
            printf("lab = "); fmprb_print(lab); printf("\n\n");
            printf("lalb = "); fmprb_print(lalb); printf("\n\n");
            abort();
        }

        fmprb_log(a, a, prec);
        if (!fmprb_overlaps(a, la))
        {
            printf("FAIL: aliasing\n\n");
            abort();
        }

        fmprb_clear(a);
        fmprb_clear(ab);
        fmprb_clear(lab);
        fmprb_clear(la);
        fmprb_clear(lb);
        fmprb_clear(lalb);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

