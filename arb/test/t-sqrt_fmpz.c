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

#include "arb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("sqrt_fmpz....");
    fflush(stdout);

    flint_randinit(state);

    /* check (sqrt(n))^2 contains n */
    for (iter = 0; iter < 10000; iter++)
    {
        arb_t a, b;
        fmpz_t n;

        arb_init(a, 1 + n_randint(state, 500));
        arb_init(b, 1 + n_randint(state, 500));

        fmpz_init(n);
        fmpz_randtest_unsigned(n, state, 1 + n_randint(state, 500));

        arb_randtest(a, state, 10);
        arb_randtest(b, state, 10);

        arb_sqrt_fmpz(a, n);
        arb_mul(b, a, a);

        if (!arb_contains_fmpz(b, n))
        {
            printf("FAIL: containment\n\n");
            printf("n = "); fmpz_print(n); printf("\n\n");
            printf("a = "); arb_debug(a); printf("\n\n");
            printf("b = "); arb_debug(b); printf("\n\n");
            abort();
        }

        arb_clear(a);
        arb_clear(b);
        fmpz_clear(n);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
