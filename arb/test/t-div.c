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

    printf("div....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        arb_t a, b, c;
        fmpq_t x, y, z;

        arb_init(a);
        arb_init(b);
        arb_init(c);

        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);

        do {
            arb_randtest(a, state, 1 + n_randint(state, 200), 10);
            arb_randtest(b, state, 1 + n_randint(state, 200), 10);
            arb_randtest(c, state, 1 + n_randint(state, 200), 10);

            arb_get_rand_fmpq(x, state, a, 1 + n_randint(state, 200));
            arb_get_rand_fmpq(y, state, b, 1 + n_randint(state, 200));
        } while (fmpq_is_zero(y));

        arb_div(c, a, b, 2 + n_randint(state, 200));
        fmpq_div(z, x, y);

        if (!arb_contains_fmpq(c, z))
        {
            printf("FAIL: containment\n\n");
            printf("a = "); arb_print(a); printf("\n\n");
            printf("x = "); fmpq_print(x); printf("\n\n");
            printf("b = "); arb_print(b); printf("\n\n");
            printf("y = "); fmpq_print(y); printf("\n\n");
            printf("c = "); arb_print(c); printf("\n\n");
            printf("z = "); fmpq_print(z); printf("\n\n");
            abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);

        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_clear(z);
    }

    /* aliasing of c and a */
    for (iter = 0; iter < 10000; iter++)
    {
        arb_t a, b;
        fmpq_t x, y, z;

        arb_init(a);
        arb_init(b);

        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);

        do {
            arb_randtest(a, state, 1 + n_randint(state, 200), 10);
            arb_randtest(b, state, 1 + n_randint(state, 200), 10);

            arb_get_rand_fmpq(x, state, a, 1 + n_randint(state, 200));
            arb_get_rand_fmpq(y, state, b, 1 + n_randint(state, 200));
        } while (fmpq_is_zero(y));

        arb_div(a, a, b, 2 + n_randint(state, 200));
        fmpq_div(z, x, y);

        if (!arb_contains_fmpq(a, z))
        {
            printf("FAIL: aliasing (c, a)\n\n");
            printf("a = "); arb_print(a); printf("\n\n");
            printf("x = "); fmpq_print(x); printf("\n\n");
            printf("b = "); arb_print(b); printf("\n\n");
            printf("y = "); fmpq_print(y); printf("\n\n");
            printf("z = "); fmpq_print(z); printf("\n\n");
            abort();
        }

        arb_clear(a);
        arb_clear(b);

        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_clear(z);
    }

    /* aliasing of c and b */
    for (iter = 0; iter < 10000; iter++)
    {
        arb_t a, b;
        fmpq_t x, y, z;

        arb_init(a);
        arb_init(b);

        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);

        do {
            arb_randtest(a, state, 1 + n_randint(state, 200), 10);
            arb_randtest(b, state, 1 + n_randint(state, 200), 10);

            arb_get_rand_fmpq(x, state, a, 1 + n_randint(state, 200));
            arb_get_rand_fmpq(y, state, b, 1 + n_randint(state, 200));
        } while (fmpq_is_zero(y));

        arb_div(b, a, b, 2 + n_randint(state, 200));
        fmpq_div(z, x, y);

        if (!arb_contains_fmpq(b, z))
        {
            printf("FAIL: aliasing (c, b)\n\n");
            printf("a = "); arb_print(a); printf("\n\n");
            printf("x = "); fmpq_print(x); printf("\n\n");
            printf("b = "); arb_print(b); printf("\n\n");
            printf("y = "); fmpq_print(y); printf("\n\n");
            printf("z = "); fmpq_print(z); printf("\n\n");
            abort();
        }

        arb_clear(a);
        arb_clear(b);

        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_clear(z);
    }

    /* test special values */
    for (iter = 0; iter < 100000; iter++)
    {
        arb_t a, b, c, d;

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arb_init(d);

        arb_randtest_special(a, state, 1 + n_randint(state, 200), 1 + n_randint(state, 100));
        arb_randtest_special(b, state, 1 + n_randint(state, 200), 1 + n_randint(state, 100));
        arb_randtest_special(c, state, 1 + n_randint(state, 200), 1 + n_randint(state, 100));

        arb_div(c, a, b, 2 + n_randint(state, 200));
        arb_mul(d, c, b, 2 + n_randint(state, 200));

        if (!arb_contains(d, a))
        {
            printf("FAIL: containment\n\n");
            printf("a = "); arb_printd(a, 15); printf("\n\n");
            printf("b = "); arb_printd(b, 15); printf("\n\n");
            printf("c = "); arb_printd(c, 15); printf("\n\n");
            printf("d = "); arb_printd(d, 15); printf("\n\n");
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
