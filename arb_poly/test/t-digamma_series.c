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

#include "arb_poly.h"

int main()
{
    slong iter;
    flint_rand_t state;

    printf("digamma_series....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        slong m, n1, n2, rbits1, rbits2, rbits3;
        arb_poly_t a, b, c, d;

        rbits1 = 2 + n_randint(state, 300);
        rbits2 = 2 + n_randint(state, 300);
        rbits3 = 2 + n_randint(state, 300);

        m = n_randint(state, 25);
        n1 = n_randint(state, 25);
        n2 = n_randint(state, 25);

        arb_poly_init(a);
        arb_poly_init(b);
        arb_poly_init(c);
        arb_poly_init(d);

        arb_poly_randtest(a, state, m, rbits1, 10);

        arb_poly_digamma_series(b, a, n1, rbits2);
        arb_poly_digamma_series(c, a, n2, rbits3);

        arb_poly_set(d, b);
        arb_poly_truncate(d, FLINT_MIN(n1, n2));
        arb_poly_truncate(c, FLINT_MIN(n1, n2));

        if (!arb_poly_overlaps(c, d))
        {
            printf("FAIL\n\n");
            printf("n1 = %ld, n2 = %ld, bits2 = %ld, bits3 = %ld\n", n1, n2, rbits2, rbits3);

            printf("a = "); arb_poly_printd(a, 50); printf("\n\n");
            printf("b = "); arb_poly_printd(b, 50); printf("\n\n");
            printf("c = "); arb_poly_printd(c, 50); printf("\n\n");

            abort();
        }

        /* check psi(a) + 1/a = psi(a+1) */
        arb_poly_inv_series(c, a, n1, rbits2);
        arb_poly_add(c, b, c, rbits2);

        arb_poly_add_si(d, a, 1, rbits2);
        arb_poly_digamma_series(d, d, n1, rbits2);

        if (!arb_poly_overlaps(c, d))
        {
            printf("FAIL (functional equation)\n\n");

            printf("a = "); arb_poly_printd(a, 15); printf("\n\n");
            printf("b = "); arb_poly_printd(b, 15); printf("\n\n");
            printf("c = "); arb_poly_printd(c, 15); printf("\n\n");
            printf("d = "); arb_poly_printd(d, 15); printf("\n\n");

            abort();
        }

        arb_poly_digamma_series(a, a, n1, rbits2);
        if (!arb_poly_overlaps(a, b))
        {
            printf("FAIL (aliasing)\n\n");
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

