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

    printf("sin....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        fmprb_t a, b;
        fmpq_t q;
        mpfr_t t;
        long prec = 2 + n_randint(state, 200);

        fmprb_init(a);
        fmprb_init(b);
        fmpq_init(q);
        mpfr_init2(t, prec + 100);

        fmprb_randtest(a, state, 1 + n_randint(state, 200), 3);
        fmprb_randtest(b, state, 1 + n_randint(state, 200), 3);
        fmprb_get_rand_fmpq(q, state, a, 1 + n_randint(state, 200));

        fmpq_get_mpfr(t, q, MPFR_RNDN);
        mpfr_sin(t, t, MPFR_RNDN);

        fmprb_sin(b, a, prec);

        if (!fmprb_contains_mpfr(b, t))
        {
            printf("FAIL: containment\n\n");
            printf("a = "); fmprb_print(a); printf("\n\n");
            printf("b = "); fmprb_print(b); printf("\n\n");
            abort();
        }

        fmprb_sin(a, a, prec);

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

    /* check large arguments */
    for (iter = 0; iter < 10000; iter++)
    {
        fmprb_t a, b, c, d;
        long prec1, prec2;

        prec1 = 2 + n_randint(state, 1000);
        prec2 = prec1 + 30;

        fmprb_init(a);
        fmprb_init(b);
        fmprb_init(c);
        fmprb_init(d);

        fmprb_randtest_precise(a, state, 1 + n_randint(state, 1000), 100);

        fmprb_sin(b, a, prec1);
        fmprb_sin(c, a, prec2);

        if (!fmprb_overlaps(b, c))
        {
            printf("FAIL: overlap\n\n");
            printf("a = "); fmprb_print(a); printf("\n\n");
            printf("b = "); fmprb_print(b); printf("\n\n");
            printf("c = "); fmprb_print(c); printf("\n\n");
            abort();
        }

        /* check sin(2a) = 2sin(a)cos(a) */
        fmprb_mul_2exp_si(c, a, 1);
        fmprb_sin(c, c, prec1);

        fmprb_cos(d, a, prec1);
        fmprb_mul(b, b, d, prec1);
        fmprb_mul_2exp_si(b, b, 1);

        if (!fmprb_overlaps(b, c))
        {
            printf("FAIL: functional equation\n\n");
            printf("a = "); fmprb_print(a); printf("\n\n");
            printf("b = "); fmprb_print(b); printf("\n\n");
            printf("c = "); fmprb_print(c); printf("\n\n");
            abort();
        }

        fmprb_clear(a);
        fmprb_clear(b);
        fmprb_clear(c);
        fmprb_clear(d);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
