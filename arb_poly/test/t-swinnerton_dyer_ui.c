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

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include "arb_poly.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("swinnerton_dyer_ui....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 50; iter++)
    {
        arb_poly_t a, b;
        arb_t x, y;
        fmpz_poly_t c;
        long i, n, prec;

        arb_poly_init(a);
        arb_poly_init(b);
        fmpz_poly_init(c);
        arb_init(x);
        arb_init(y);

        n = n_randint(state, 10);
        arb_poly_swinnerton_dyer_ui(a, n, 0);
        prec = 2 + n_randint(state, 10000);

        if (!arb_poly_get_unique_fmpz_poly(c, a))
        {
            printf("FAIL (uniqueness)\n\n");
            abort();
        }

        arb_poly_set_fmpz_poly(b, c, prec);

        arb_zero(x);
        for (i = 0; i < n; i++)
        {
            arb_sqrt_ui(y, n_nth_prime(i + 1), prec);
            if (n_randint(state, 2))
                arb_add(x, x, y, prec);
            else
                arb_sub(x, x, y, prec);
        }

        arb_poly_evaluate(y, b, x, prec);

        if (!arb_contains_zero(y))
        {
            printf("FAIL (containment)\n\n");
            abort();
        }

        arb_poly_swinnerton_dyer_ui(b, n, 2 + n_randint(state, 1000));

        if (!arb_poly_overlaps(a, b))
        {
            printf("FAIL (overlap)\n\n");
            abort();
        }

        arb_poly_clear(a);
        arb_poly_clear(b);
        fmpz_poly_clear(c);
        arb_clear(x);
        arb_clear(y);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

