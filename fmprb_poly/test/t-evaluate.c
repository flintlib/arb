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

#include "fmprb_poly.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("evaluate....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        long qbits1, qbits2, rbits1, rbits2, rbits3;
        fmpq_poly_t F;
        fmpq_t X, Y;
        fmprb_poly_t f;
        fmprb_t x, y;

        qbits1 = 2 + n_randint(state, 200);
        qbits2 = 2 + n_randint(state, 200);
        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);
        rbits3 = 2 + n_randint(state, 200);

        fmpq_poly_init(F);
        fmpq_init(X);
        fmpq_init(Y);

        fmprb_poly_init(f);
        fmprb_init(x);
        fmprb_init(y);

        fmpq_poly_randtest(F, state, 1 + n_randint(state, 20), qbits1);
        fmpq_randtest(X, state, qbits2);
        fmpq_poly_evaluate_fmpq(Y, F, X);

        fmprb_poly_set_fmpq_poly(f, F, rbits1);
        fmprb_set_fmpq(x, X, rbits2);
        fmprb_poly_evaluate(y, f, x, rbits3);

        if (!fmprb_contains_fmpq(y, Y))
        {
            printf("FAIL\n\n");

            printf("F = "); fmpq_poly_print(F); printf("\n\n");
            printf("X = "); fmpq_print(X); printf("\n\n");
            printf("Y = "); fmpq_print(Y); printf("\n\n");

            printf("f = "); fmprb_poly_printd(f, 15); printf("\n\n");
            printf("x = "); fmprb_printd(x, 15); printf("\n\n");
            printf("y = "); fmprb_printd(y, 15); printf("\n\n");

            abort();
        }

        /* aliasing */
        fmprb_poly_evaluate(x, f, x, rbits3);
        if (!fmprb_contains_fmpq(x, Y))
        {
            printf("FAIL (aliasing)\n\n");
            abort();
        }

        fmpq_poly_clear(F);
        fmpq_clear(X);
        fmpq_clear(Y);

        fmprb_poly_clear(f);
        fmprb_clear(x);
        fmprb_clear(y);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
