/*=============================================================================

    This file is part of ARB.

    ARB is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    ARB is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ARB; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "arb_poly.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("set_fmpq_poly....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        arb_poly_t a;
        fmpq_poly_t x;

        arb_poly_init(a, 1 + n_randint(state, 200));
        fmpq_poly_init(x);

        arb_poly_randtest(a, state, n_randint(state, 20), 10);
        fmpq_poly_randtest(x, state, n_randint(state, 20), 1 + n_randint(state, 200));
        arb_poly_set_fmpq_poly(a, x);

        if (!arb_poly_contains_fmpq_poly(a, x))
        {
            printf("FAIL\n\n");
            printf("a = "); arb_poly_debug(a); printf("\n\n");
            printf("x = "); fmpq_poly_print(x); printf("\n\n");
            abort();
        }

        arb_poly_clear(a);
        fmpq_poly_clear(x);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
