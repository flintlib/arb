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

#include "fmprb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("acos....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
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
        mpfr_acos(t, t, MPFR_RNDN);

        fmprb_acos(b, a, prec);

        if (!fmprb_contains_mpfr(b, t))
        {
            printf("FAIL: containment, prec = %ld\n\n", prec);
            printf("a = "); fmprb_printd(a, 100); printf("\n\n");
            printf("b = "); fmprb_printd(b, 100); printf("\n\n");
            abort();
        }

        fmprb_acos(a, a, prec);

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
        fmprb_t a, b, c;
        long prec1, prec2;

        prec1 = 2 + n_randint(state, 1000);
        prec2 = prec1 + 30;

        fmprb_init(a);
        fmprb_init(b);
        fmprb_init(c);

        fmprb_randtest_precise(a, state, 1 + n_randint(state, 1000), 100);

        fmprb_acos(b, a, prec1);
        fmprb_acos(c, a, prec2);

        if (!fmprb_overlaps(b, c))
        {
            printf("FAIL: overlap\n\n");
            printf("a = "); fmprb_print(a); printf("\n\n");
            printf("b = "); fmprb_print(b); printf("\n\n");
            printf("c = "); fmprb_print(c); printf("\n\n");
            abort();
        }

        /* check sin(asin(x)) = x */
        fmprb_cos(c, b, prec1);

        if (!fmprb_contains(c, a))
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
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
