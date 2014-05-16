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

    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include "acb_poly.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("get_coeff_ptr....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 1000; i++)
    {
        acb_poly_t A;
        acb_t a;
        long n = n_randint(state, 100);

        acb_poly_init(A);
        acb_poly_randtest(A, state, n_randint(state, 100), 100, 10);
        acb_init(a);

        acb_poly_get_coeff_acb(a, A, n);

        result = n < acb_poly_length(A) ? 
                     acb_equal(a, acb_poly_get_coeff_ptr(A, n)) : 
                     acb_poly_get_coeff_ptr(A, n) == NULL;
        if (!result)
        {
            printf("FAIL:\n");
            printf("A = "), acb_poly_printd(A, 10), printf("\n\n");
            printf("a = "), acb_print(a), printf("\n\n");
            printf("n = %ld\n\n", n);
            abort();
        }

        acb_poly_clear(A);
        acb_clear(a);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return 0;
}

