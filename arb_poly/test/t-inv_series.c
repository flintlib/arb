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

#include "arb_poly.h"


int main()
{
    long iter;
    flint_rand_t state;

    printf("inv_series....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        long m, n, qbits, rbits1, rbits2;
        fmpq_poly_t A, B;
        arb_poly_t a, b;

        qbits = 2 + n_randint(state, 200);
        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);

        m = 1 + n_randint(state, 20);
        n = 1 + n_randint(state, 20);

        fmpq_poly_init(A);
        fmpq_poly_init(B);

        arb_poly_init(a);
        arb_poly_init(b);

        do {
            fmpq_poly_randtest_not_zero(A, state, m, qbits);
        } while (A->coeffs[0] == 0);

        arb_poly_randtest(b, state, 1 + n_randint(state, 20), rbits1, 5);

        fmpq_poly_inv_series(B, A, n);
        arb_poly_set_fmpq_poly(a, A, rbits1);
        arb_poly_inv_series(b, a, n, rbits2);

        if (!arb_poly_contains_fmpq_poly(b, B))
        {
            printf("FAIL\n\n");
            printf("bits2 = %ld\n", rbits2);

            printf("A = "); fmpq_poly_print(A); printf("\n\n");
            printf("B = "); fmpq_poly_print(B); printf("\n\n");

            printf("a = "); arb_poly_printd(a, 15); printf("\n\n");
            printf("b = "); arb_poly_printd(b, 15); printf("\n\n");

            abort();
        }

        arb_poly_inv_series(a, a, n, rbits2);
        if (!arb_poly_equal(a, b))
        {
            printf("FAIL (aliasing)\n\n");
            abort();
        }

        fmpq_poly_clear(A);
        fmpq_poly_clear(B);

        arb_poly_clear(a);
        arb_poly_clear(b);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
