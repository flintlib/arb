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

    printf("sin_cos....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        arb_t a, b, c;
        fmpq_t q;
        mpfr_t t, u;
        long prec0, prec;

        prec0 = 400;
        if (iter % 100 == 0)
            prec0 = 8000;

        prec = 2 + n_randint(state, prec0);

        arb_init(a);
        arb_init(b);
        arb_init(c);
        fmpq_init(q);
        mpfr_init2(t, prec + 100);
        mpfr_init2(u, prec + 100);

        arb_randtest(a, state, 1 + n_randint(state, prec0), 6);
        arb_randtest(b, state, 1 + n_randint(state, prec0), 6);
        arb_randtest(c, state, 1 + n_randint(state, prec0), 6);
        arb_get_rand_fmpq(q, state, a, 1 + n_randint(state, prec0));

        fmpq_get_mpfr(t, q, MPFR_RNDN);
        mpfr_sin_cos(t, u, t, MPFR_RNDN);

        arb_sin_cos(b, c, a, prec);

        if (!arb_contains_mpfr(b, t))
        {
            printf("FAIL: containment (sin)\n\n");
            printf("a = "); arb_print(a); printf("\n\n");
            printf("b = "); arb_print(b); printf("\n\n");
            abort();
        }

        if (!arb_contains_mpfr(c, u))
        {
            printf("FAIL: containment (cos)\n\n");
            printf("a = "); arb_print(a); printf("\n\n");
            printf("c = "); arb_print(c); printf("\n\n");
            abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        fmpq_clear(q);
        mpfr_clear(t);
        mpfr_clear(u);
    }

    /* check large arguments */
    for (iter = 0; iter < 1000000; iter++)
    {
        arb_t a, b, c, d, e;
        long prec0, prec1, prec2, prec3;

        if (iter % 10 == 0)
            prec0 = 10000;
        else
            prec0 = 1000;

        prec1 = 2 + n_randint(state, prec0);
        prec2 = 2 + n_randint(state, prec0);

        if (iter % 10 == 0)
            prec3 = 50000;
        else
            prec3 = 100;

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arb_init(d);
        arb_init(e);

        arb_randtest_special(a, state, 1 + n_randint(state, prec0), prec3);
        arb_randtest_special(b, state, 1 + n_randint(state, prec0), 100);
        arb_randtest_special(c, state, 1 + n_randint(state, prec0), 100);
        arb_randtest_special(d, state, 1 + n_randint(state, prec0), 100);
        arb_randtest_special(e, state, 1 + n_randint(state, prec0), 100);

        arb_sin_cos(b, c, a, prec1);
        arb_sin_cos(d, e, a, prec2);

        if (!arb_overlaps(b, d) || !arb_overlaps(c, e))
        {
            printf("FAIL: overlap\n\n");
            printf("a = "); arb_print(a); printf("\n\n");
            printf("b = "); arb_print(b); printf("\n\n");
            printf("c = "); arb_print(c); printf("\n\n");
            printf("d = "); arb_print(d); printf("\n\n");
            printf("e = "); arb_print(e); printf("\n\n");
            abort();
        }

        /* check sin(a)^2 + cos(a)^2 = 1 */
        arb_mul(d, b, b, prec1);
        arb_mul(e, c, c, prec1);
        arb_add(d, d, e, prec1);
        arb_sub_ui(d, d, 1, prec1);

        if (!arb_contains_zero(d))
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
        arb_clear(d);
        arb_clear(e);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
