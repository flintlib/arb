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

#include "fmpcb_poly.h"


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
        fmpcb_poly_t a, b;

        qbits = 2 + n_randint(state, 200);
        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);

        m = 1 + n_randint(state, 20);
        n = 1 + n_randint(state, 20);

        fmpq_poly_init(A);
        fmpq_poly_init(B);

        fmpcb_poly_init(a);
        fmpcb_poly_init(b);

        do {
            fmpq_poly_randtest_not_zero(A, state, m, qbits);
        } while (A->coeffs[0] == 0);

        fmpq_poly_inv_series(B, A, n);
        fmpcb_poly_set_fmpq_poly(a, A, rbits1);
        fmpcb_poly_inv_series(b, a, n, rbits2);

        if (!fmpcb_poly_contains_fmpq_poly(b, B))
        {
            printf("FAIL\n\n");
            printf("bits2 = %ld\n", rbits2);

            printf("A = "); fmpq_poly_print(A); printf("\n\n");
            printf("B = "); fmpq_poly_print(B); printf("\n\n");

            printf("a = "); fmpcb_poly_printd(a, 15); printf("\n\n");
            printf("b = "); fmpcb_poly_printd(b, 15); printf("\n\n");

            abort();
        }

        fmpcb_poly_inv_series(a, a, n, rbits2);
        if (!fmpcb_poly_equal(a, b))
        {
            printf("FAIL (aliasing)\n\n");
            abort();
        }

        fmpq_poly_clear(A);
        fmpq_poly_clear(B);

        fmpcb_poly_clear(a);
        fmpcb_poly_clear(b);
    }

    /* check f * f^-1 = 1 */
    for (iter = 0; iter < 10000; iter++)
    {
        long bits, trunc;
        fmpcb_poly_t a, b, ab;
        fmpq_poly_t id;

        bits = 2 + n_randint(state, 200);
        trunc = 1 + n_randint(state, 10);

        fmpcb_poly_init(a);
        fmpcb_poly_init(b);
        fmpcb_poly_init(ab);
        fmpq_poly_init(id);

        do {
            fmpcb_poly_randtest(a, state, 1 + n_randint(state, 10), bits, 5);
        } while (a->length == 0 || fmpcb_contains_zero(a->coeffs + 0));

        fmpcb_poly_inv_series(b, a, trunc, bits);
        fmpcb_poly_mullow(ab, a, b, trunc, bits);

        fmpq_poly_set_ui(id, 1);

        if (!fmpcb_poly_contains_fmpq_poly(ab, id))
        {
            printf("FAIL\n\n");
            printf("bits = %ld\n", bits);
            printf("trunc = %ld\n", trunc);

            printf("a = "); fmpcb_poly_printd(a, 15); printf("\n\n");
            printf("b = "); fmpcb_poly_printd(b, 15); printf("\n\n");
            printf("ab = "); fmpcb_poly_printd(ab, 15); printf("\n\n");

            abort();
        }

        fmpcb_poly_clear(a);
        fmpcb_poly_clear(b);
        fmpcb_poly_clear(ab);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
