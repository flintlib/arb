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

    printf("pow_ui....");
    fflush(stdout);

    flint_randinit(state);

    /* compare with fmpz_poly */
    for (iter = 0; iter < 10000; iter++)
    {
        long zbits1, rbits1, rbits2;
        ulong e;
        fmpz_poly_t A, B;
        arb_poly_t a, b;

        zbits1 = 2 + n_randint(state, 100);
        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);
        e = n_randint(state, 30);

        fmpz_poly_init(A);
        fmpz_poly_init(B);

        arb_poly_init(a);
        arb_poly_init(b);

        fmpz_poly_randtest(A, state, 1 + n_randint(state, 8), zbits1);
        fmpz_poly_pow(B, A, e);

        arb_poly_set_fmpz_poly(a, A, rbits1);
        arb_poly_pow_ui(b, a, e, rbits2);

        if (!arb_poly_contains_fmpz_poly(b, B))
        {
            printf("FAIL\n\n");
            printf("bits2 = %ld\n", rbits2);
            printf("e = %lu\n", e);

            printf("A = "); fmpz_poly_print(A); printf("\n\n");
            printf("B = "); fmpz_poly_print(B); printf("\n\n");

            printf("a = "); arb_poly_printd(a, 15); printf("\n\n");
            printf("b = "); arb_poly_printd(b, 15); printf("\n\n");

            abort();
        }

        arb_poly_pow_ui(a, a, e, rbits2);
        if (!arb_poly_equal(a, b))
        {
            printf("FAIL (aliasing)\n\n");
            abort();
        }

        fmpz_poly_clear(A);
        fmpz_poly_clear(B);

        arb_poly_clear(a);
        arb_poly_clear(b);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
