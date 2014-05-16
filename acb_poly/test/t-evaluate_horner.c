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

#include "acb_poly.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("evaluate_horner....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        long qbits1, qbits2, rbits1, rbits2, rbits3;
        fmpq_poly_t F;
        fmpq_t X, Y;
        acb_poly_t f;
        acb_t x, y;

        qbits1 = 2 + n_randint(state, 200);
        qbits2 = 2 + n_randint(state, 200);
        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);
        rbits3 = 2 + n_randint(state, 200);

        fmpq_poly_init(F);
        fmpq_init(X);
        fmpq_init(Y);

        acb_poly_init(f);
        acb_init(x);
        acb_init(y);

        fmpq_poly_randtest(F, state, 1 + n_randint(state, 20), qbits1);
        fmpq_randtest(X, state, qbits2);
        fmpq_poly_evaluate_fmpq(Y, F, X);

        acb_poly_set_fmpq_poly(f, F, rbits1);
        acb_set_fmpq(x, X, rbits2);
        acb_poly_evaluate_horner(y, f, x, rbits3);

        if (!acb_contains_fmpq(y, Y))
        {
            printf("FAIL\n\n");

            printf("F = "); fmpq_poly_print(F); printf("\n\n");
            printf("X = "); fmpq_print(X); printf("\n\n");
            printf("Y = "); fmpq_print(Y); printf("\n\n");

            printf("f = "); acb_poly_printd(f, 15); printf("\n\n");
            printf("x = "); acb_printd(x, 15); printf("\n\n");
            printf("y = "); acb_printd(y, 15); printf("\n\n");

            abort();
        }

        /* aliasing */
        acb_poly_evaluate_horner(x, f, x, rbits3);
        if (!acb_contains_fmpq(x, Y))
        {
            printf("FAIL (aliasing)\n\n");
            abort();
        }

        fmpq_poly_clear(F);
        fmpq_clear(X);
        fmpq_clear(Y);

        acb_poly_clear(f);
        acb_clear(x);
        acb_clear(y);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
