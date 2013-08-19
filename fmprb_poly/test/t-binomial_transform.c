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

#include "fmprb_poly.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("binomial_transform....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 5000; iter++)
    {
        fmprb_poly_t a, b, c, d;
        long j, n, prec;

        fmprb_poly_init(a);
        fmprb_poly_init(b);
        fmprb_poly_init(c);
        fmprb_poly_init(d);

        n = n_randint(state, 20);
        prec = 2 + n_randint(state, 200);

        fmprb_poly_randtest(a, state, n, prec, 10);
        fmprb_poly_randtest(b, state, n, prec, 10);
        fmprb_poly_randtest(c, state, n, prec, 10);

        /* check self-inversion property */
        fmprb_poly_binomial_transform(b, a, n, prec);
        fmprb_poly_binomial_transform(c, b, n, prec);

        fmprb_poly_set(d, a);
        fmprb_poly_truncate(d, n);

        if (!fmprb_poly_contains(c, d))
        {
            printf("FAIL (containment)\n\n");
            printf("n = %ld, prec = %ld\n\n", n, prec);

            printf("a: "); fmprb_poly_printd(a, 15); printf("\n\n");
            printf("b: "); fmprb_poly_printd(b, 15); printf("\n\n");
            printf("c: "); fmprb_poly_printd(c, 15); printf("\n\n");
            printf("d: "); fmprb_poly_printd(d, 15); printf("\n\n");

            abort();
        }

        fmprb_poly_set(d, a);
        fmprb_poly_binomial_transform(d, d, n, prec);
        if (!fmprb_poly_equal(d, b))
        {
            printf("FAIL (aliasing)\n\n");

            printf("a: "); fmprb_poly_printd(a, 15); printf("\n\n");
            printf("b: "); fmprb_poly_printd(b, 15); printf("\n\n");
            printf("d: "); fmprb_poly_printd(d, 15); printf("\n\n");

            abort();
        }

        /* compare with power series operations */
        fmprb_poly_zero(d);
        for (j = 1; j < n; j++)
            fmprb_poly_set_coeff_si(d, j, -1);
        fmprb_poly_compose_series(c, a, d, n, prec);
        for (j = 0; j < n; j++)
            fmprb_poly_set_coeff_si(d, j, 1);
        fmprb_poly_mullow(c, c, d, n, prec);

        if (!fmprb_poly_overlaps(b, c))
        {
            printf("FAIL (power series)\n\n");

            printf("a: "); fmprb_poly_printd(a, 15); printf("\n\n");
            printf("b: "); fmprb_poly_printd(b, 15); printf("\n\n");
            printf("c: "); fmprb_poly_printd(c, 15); printf("\n\n");

            abort();
        }

        fmprb_poly_clear(a);
        fmprb_poly_clear(b);
        fmprb_poly_clear(c);
        fmprb_poly_clear(d);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

