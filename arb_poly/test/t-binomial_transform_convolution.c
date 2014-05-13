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

    printf("binomial_transform_convolution....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 5000; iter++)
    {
        arb_poly_t a, b, c, d;
        long j, n, prec;

        arb_poly_init(a);
        arb_poly_init(b);
        arb_poly_init(c);
        arb_poly_init(d);

        n = n_randint(state, 20);
        prec = 2 + n_randint(state, 200);

        arb_poly_randtest(a, state, n, prec, 10);
        arb_poly_randtest(b, state, n, prec, 10);
        arb_poly_randtest(c, state, n, prec, 10);

        /* check self-inversion property */
        arb_poly_binomial_transform_convolution(b, a, n, prec);
        arb_poly_binomial_transform_convolution(c, b, n, prec);

        arb_poly_set(d, a);
        arb_poly_truncate(d, n);

        if (!arb_poly_contains(c, d))
        {
            printf("FAIL (containment)\n\n");
            printf("n = %ld, prec = %ld\n\n", n, prec);

            printf("a: "); arb_poly_printd(a, 15); printf("\n\n");
            printf("b: "); arb_poly_printd(b, 15); printf("\n\n");
            printf("c: "); arb_poly_printd(c, 15); printf("\n\n");
            printf("d: "); arb_poly_printd(d, 15); printf("\n\n");

            abort();
        }

        arb_poly_set(d, a);
        arb_poly_binomial_transform_convolution(d, d, n, prec);
        if (!arb_poly_equal(d, b))
        {
            printf("FAIL (aliasing)\n\n");

            printf("a: "); arb_poly_printd(a, 15); printf("\n\n");
            printf("b: "); arb_poly_printd(b, 15); printf("\n\n");
            printf("d: "); arb_poly_printd(d, 15); printf("\n\n");

            abort();
        }

        /* compare with power series operations */
        arb_poly_zero(d);
        for (j = 1; j < n; j++)
            arb_poly_set_coeff_si(d, j, -1);
        arb_poly_compose_series(c, a, d, n, prec);
        for (j = 0; j < n; j++)
            arb_poly_set_coeff_si(d, j, 1);
        arb_poly_mullow(c, c, d, n, prec);

        if (!arb_poly_overlaps(b, c))
        {
            printf("FAIL (power series)\n\n");

            printf("a: "); arb_poly_printd(a, 15); printf("\n\n");
            printf("b: "); arb_poly_printd(b, 15); printf("\n\n");
            printf("c: "); arb_poly_printd(c, 15); printf("\n\n");

            abort();
        }

        arb_poly_clear(a);
        arb_poly_clear(b);
        arb_poly_clear(c);
        arb_poly_clear(d);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

