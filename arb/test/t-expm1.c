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

    Copyright (C) 2012, 2013 Fredrik Johansson

******************************************************************************/

#include "arb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("expm1....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        arb_t a, b;
        fmpq_t q;
        mpfr_t t;
        long prec0, prec;

        prec0 = 400;
        if (iter % 100 == 0)
            prec0 = 10000;

        prec = 2 + n_randint(state, prec0);

        arb_init(a);
        arb_init(b);
        fmpq_init(q);
        mpfr_init2(t, prec + 100);

        arb_randtest(a, state, 1 + n_randint(state, prec0), 4);
        arb_randtest(b, state, 1 + n_randint(state, prec0), 4);
        arb_get_rand_fmpq(q, state, a, 1 + n_randint(state, prec0));

        fmpq_get_mpfr(t, q, MPFR_RNDN);
        mpfr_expm1(t, t, MPFR_RNDN);

        arb_expm1(b, a, prec);

        if (!arb_contains_mpfr(b, t))
        {
            printf("FAIL: containment\n\n");
            printf("iter = %ld, prec = %ld\n\n", iter, prec);
            printf("a = "); arb_print(a); printf("\n\n");
            printf("b = "); arb_print(b); printf("\n\n");
            abort();
        }

        arb_expm1(a, a, prec);

        if (!arb_equal(a, b))
        {
            printf("FAIL: aliasing\n\n");
            printf("iter = %ld, prec = %ld\n\n", iter, prec);
            printf("a = "); arb_print(a); printf("\n\n");
            printf("b = "); arb_print(b); printf("\n\n");
            abort();
        }

        arb_clear(a);
        arb_clear(b);
        fmpq_clear(q);
        mpfr_clear(t);
    }

    /* check large arguments */
    for (iter = 0; iter < 100000; iter++)
    {
        arb_t a, b, c, d;
        long prec0, prec1, prec2;

        if (iter % 10 == 0)
            prec0 = 10000;
        else
            prec0 = 1000;

        prec1 = 2 + n_randint(state, prec0);
        prec2 = 2 + n_randint(state, prec0);

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arb_init(d);

        arb_randtest_special(a, state, 1 + n_randint(state, prec0), 100);
        arb_randtest_special(b, state, 1 + n_randint(state, prec0), 100);

        arb_expm1(b, a, prec1);
        arb_expm1(c, a, prec2);

        if (!arb_overlaps(b, c))
        {
            printf("FAIL: overlap\n\n");
            printf("a = "); arb_print(a); printf("\n\n");
            printf("b = "); arb_print(b); printf("\n\n");
            printf("c = "); arb_print(c); printf("\n\n");
            abort();
        }

        arb_randtest_special(b, state, 1 + n_randint(state, prec0), 100);

        /* compare with exp */
        arb_expm1(c, a, prec1);

        arb_exp(d, a, prec1);
        arb_sub_ui(d, d, 1, prec1);

        if (!arb_overlaps(c, d))
        {
            printf("FAIL: comparison with exp\n\n");
            printf("a = "); arb_print(a); printf("\n\n");
            printf("b = "); arb_print(b); printf("\n\n");
            printf("c = "); arb_print(c); printf("\n\n");
            printf("d = "); arb_print(d); printf("\n\n");
            abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        arb_clear(d);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
