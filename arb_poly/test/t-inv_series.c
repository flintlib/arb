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

static int
is_invertible(arb_poly_t poly)
{
    if (poly->length == 0)
        return 0;

    return !(fmpz_cmpabs(arb_poly_radref(poly), arb_poly_coeffs(poly)) >= 0);
}

int main()
{
    long iter;
    flint_rand_t state;

    printf("inv_series....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        arb_poly_t a, b;
        fmpq_poly_t x, y;
        long n;

        arb_poly_init(a, 1 + n_randint(state, 200));
        arb_poly_init(b, 1 + n_randint(state, 200));

        fmpq_poly_init(x);
        fmpq_poly_init(y);

        do {
            arb_poly_randtest(a, state, n_randint(state, 20), 9);
        } while (!is_invertible(a));

        arb_poly_randtest(b, state, n_randint(state, 20), 9);

        arb_poly_get_rand_fmpq_poly(x, state, a);
        arb_poly_get_rand_fmpq_poly(y, state, b);

        n = 1 + n_randint(state, 20);

        arb_poly_inv_series(b, a, n);
        fmpq_poly_inv_series(y, x, n);

        if (!arb_poly_contains_fmpq_poly(b, y))
        {
            printf("FAIL: containment\n\n");
            printf("a = "); arb_poly_debug(a); printf("\n\n");
            printf("x = "); fmpq_poly_print(x); printf("\n\n");
            printf("b = "); arb_poly_debug(b); printf("\n\n");
            printf("y = "); fmpq_poly_print(y); printf("\n\n");
            abort();
        }

        arb_poly_clear(a);
        arb_poly_clear(b);

        fmpq_poly_clear(x);
        fmpq_poly_clear(y);
    }

    /* aliasing */
    for (iter = 0; iter < 10000; iter++)
    {
        arb_poly_t a;
        fmpq_poly_t x;
        long n;

        arb_poly_init(a, 1 + n_randint(state, 200));
        fmpq_poly_init(x);

        do {
            arb_poly_randtest(a, state, n_randint(state, 20), 9);
        } while (!is_invertible(a));

        arb_poly_get_rand_fmpq_poly(x, state, a);

        n = 1 + n_randint(state, 20);

        arb_poly_inv_series(a, a, n);
        fmpq_poly_inv_series(x, x, n);

        if (!arb_poly_contains_fmpq_poly(a, x))
        {
            printf("FAIL: aliasing\n\n");
            printf("a = "); arb_poly_debug(a); printf("\n\n");
            printf("x = "); fmpq_poly_print(x); printf("\n\n");
            abort();
        }

        arb_poly_clear(a);
        fmpq_poly_clear(x);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
