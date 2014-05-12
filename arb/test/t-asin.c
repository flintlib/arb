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

    printf("asin....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        arb_t a, b;
        fmpq_t q;
        mpfr_t t;
        long prec = 2 + n_randint(state, 200);

        arb_init(a);
        arb_init(b);
        fmpq_init(q);
        mpfr_init2(t, prec + 100);

        arb_randtest(a, state, 1 + n_randint(state, 200), 3);
        arb_randtest(b, state, 1 + n_randint(state, 200), 3);
        arb_get_rand_fmpq(q, state, a, 1 + n_randint(state, 200));

        fmpq_get_mpfr(t, q, MPFR_RNDN);
        mpfr_asin(t, t, MPFR_RNDN);

        arb_asin(b, a, prec);

        if (!arb_contains_mpfr(b, t))
        {
            printf("FAIL: containment\n\n");
            printf("a = "); arb_print(a); printf("\n\n");
            printf("b = "); arb_print(b); printf("\n\n");
            abort();
        }

        arb_asin(a, a, prec);

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

    /* check large arguments */
    for (iter = 0; iter < 10000; iter++)
    {
        arb_t a, b, c;
        long prec1, prec2;

        prec1 = 2 + n_randint(state, 1000);
        prec2 = prec1 + 30;

        arb_init(a);
        arb_init(b);
        arb_init(c);

        arb_randtest_precise(a, state, 1 + n_randint(state, 1000), 100);

        arb_asin(b, a, prec1);
        arb_asin(c, a, prec2);

        if (!arb_overlaps(b, c))
        {
            printf("FAIL: overlap\n\n");
            printf("a = "); arb_print(a); printf("\n\n");
            printf("b = "); arb_print(b); printf("\n\n");
            printf("c = "); arb_print(c); printf("\n\n");
            abort();
        }

        /* check sin(asin(x)) = x */
        arb_sin(c, b, prec1);

        if (!arb_contains(c, a))
        {
            printf("FAIL: functional equation\n\n");
            printf("a = "); arb_print(a); printf("\n\n");
            printf("b = "); arb_print(b); printf("\n\n");
            printf("c = "); arb_print(c); printf("\n\n");
            abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
