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

    Copyright (C) 2015 Fredrik Johansson

******************************************************************************/

#include "acb_poly.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("get_unique_fmpz_poly....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        long prec, c;
        fmpz_poly_t A, B, C;
        acb_poly_t a, b;

        fmpz_poly_init(A);
        fmpz_poly_init(B);
        fmpz_poly_init(C);
        acb_poly_init(a);
        acb_poly_init(b);

        fmpz_poly_randtest(A, state, 1 + n_randint(state, 10), 1 + n_randint(state, 1000));
        fmpz_poly_randtest(B, state, 1 + n_randint(state, 10), 1 + n_randint(state, 1000));
        fmpz_poly_randtest(C, state, 1 + n_randint(state, 10), 1 + n_randint(state, 1000));
        c = 1 + n_randint(state, 1000);

        prec = 2 + n_randint(state, 100);

        for ( ; ; )
        {
            acb_poly_set_fmpz_poly(a, A, prec);
            acb_poly_set2_fmpz_poly(b, B, C, prec);
            acb_poly_scalar_mul_2exp_si(b, b, -c);
            acb_poly_add(a, a, b, prec);
            acb_poly_sub(a, a, b, prec);

            if (acb_poly_get_unique_fmpz_poly(B, a))
            {
                if (!fmpz_poly_equal(A, B))
                {
                    printf("FAIL\n\n");
                    printf("A = "); fmpz_poly_print(A); printf("\n\n");
                    printf("B = "); fmpz_poly_print(B); printf("\n\n");
                    printf("a = "); acb_poly_printd(a, 15); printf("\n\n");
                    abort();
                }

                break;
            }
            else
            {
                prec *= 2;
            }
        }

        fmpz_poly_clear(A);
        fmpz_poly_clear(B);
        fmpz_poly_clear(C);
        acb_poly_clear(a);
        acb_poly_clear(b);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

