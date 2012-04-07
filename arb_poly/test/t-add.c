/*=============================================================================

    This file is part of ARB.

    ARB is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    ARB is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ARB; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "arb_poly.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("add....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        arb_poly_t a, b, c;
        fmpq_poly_t x, y, z;

        arb_poly_init(a, 1 + n_randint(state, 200));
        arb_poly_init(b, 1 + n_randint(state, 200));
        arb_poly_init(c, 1 + n_randint(state, 200));

        fmpq_poly_init(x);
        fmpq_poly_init(y);
        fmpq_poly_init(z);

        arb_poly_randtest(a, state, n_randint(state, 20), 10);
        arb_poly_randtest(b, state, n_randint(state, 20), 10);
        arb_poly_randtest(c, state, n_randint(state, 20), 10);

        arb_poly_get_rand_fmpq_poly(x, state, a);
        arb_poly_get_rand_fmpq_poly(y, state, b);

        arb_poly_add(c, a, b);
        fmpq_poly_add(z, x, y);

        if (!arb_poly_contains_fmpq_poly(c, z))
        {
            printf("FAIL: containment\n\n");
            printf("a = "); arb_poly_debug(a); printf("\n\n");
            printf("x = "); fmpq_poly_print(x); printf("\n\n");
            printf("b = "); arb_poly_debug(b); printf("\n\n");
            printf("y = "); fmpq_poly_print(y); printf("\n\n");
            printf("c = "); arb_poly_debug(c); printf("\n\n");
            printf("z = "); fmpq_poly_print(z); printf("\n\n");
            abort();
        }

        arb_poly_clear(a);
        arb_poly_clear(b);
        arb_poly_clear(c);

        fmpq_poly_clear(x);
        fmpq_poly_clear(y);
        fmpq_poly_clear(z);
    }

    /* aliasing of c and a */
    for (iter = 0; iter < 10000; iter++)
    {
        arb_poly_t a, b;
        fmpq_poly_t x, y, z;

        arb_poly_init(a, 1 + n_randint(state, 200));
        arb_poly_init(b, 1 + n_randint(state, 200));

        fmpq_poly_init(x);
        fmpq_poly_init(y);
        fmpq_poly_init(z);

        arb_poly_randtest(a, state, n_randint(state, 20), 10);
        arb_poly_randtest(b, state, n_randint(state, 20), 10);

        arb_poly_get_rand_fmpq_poly(x, state, a);
        arb_poly_get_rand_fmpq_poly(y, state, b);

        arb_poly_add(a, a, b);
        fmpq_poly_add(z, x, y);

        if (!arb_poly_contains_fmpq_poly(a, z))
        {
            printf("FAIL: aliasing (c, a)\n\n");
            printf("a = "); arb_poly_debug(a); printf("\n\n");
            printf("x = "); fmpq_poly_print(x); printf("\n\n");
            printf("b = "); arb_poly_debug(b); printf("\n\n");
            printf("y = "); fmpq_poly_print(y); printf("\n\n");
            printf("z = "); fmpq_poly_print(z); printf("\n\n");
            abort();
        }

        arb_poly_clear(a);
        arb_poly_clear(b);

        fmpq_poly_clear(x);
        fmpq_poly_clear(y);
        fmpq_poly_clear(z);
    }

    /* aliasing of c and b */
    for (iter = 0; iter < 10000; iter++)
    {
        arb_poly_t a, b;
        fmpq_poly_t x, y, z;

        arb_poly_init(a, 1 + n_randint(state, 200));
        arb_poly_init(b, 1 + n_randint(state, 200));

        fmpq_poly_init(x);
        fmpq_poly_init(y);
        fmpq_poly_init(z);

        arb_poly_randtest(a, state, n_randint(state, 20), 10);
        arb_poly_randtest(b, state, n_randint(state, 20), 10);

        arb_poly_get_rand_fmpq_poly(x, state, a);
        arb_poly_get_rand_fmpq_poly(y, state, b);

        arb_poly_add(b, a, b);
        fmpq_poly_add(z, x, y);

        if (!arb_poly_contains_fmpq_poly(b, z))
        {
            printf("FAIL: aliasing (c, b)\n\n");
            printf("a = "); arb_poly_debug(a); printf("\n\n");
            printf("x = "); fmpq_poly_print(x); printf("\n\n");
            printf("b = "); arb_poly_debug(b); printf("\n\n");
            printf("y = "); fmpq_poly_print(y); printf("\n\n");
            printf("z = "); fmpq_poly_print(z); printf("\n\n");
            abort();
        }

        arb_poly_clear(a);
        arb_poly_clear(b);

        fmpq_poly_clear(x);
        fmpq_poly_clear(y);
        fmpq_poly_clear(z);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
