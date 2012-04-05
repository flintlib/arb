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

    printf("sqrt_ui....");
    fflush(stdout);

    flint_randinit(state);

    /* check (sqrt(a))^2 contains a */
    for (iter = 0; iter < 10000; iter++)
    {
        arb_t a, b;
        ulong n;

        arb_init(a, 1 + n_randint(state, 500));
        arb_init(b, 1 + n_randint(state, 500));

        n = n_randtest(state);

        arb_randtest(a, state, 10);
        arb_randtest(b, state, 10);

        arb_sqrt_ui(a, n);
        arb_mul(b, a, a);

        if (!arb_contains_ui(b, n))
        {
            printf("FAIL: containment\n\n");
            printf("n = %lu\n\n", n);
            printf("a = "); arb_debug(a); printf("\n\n");
            printf("b = "); arb_debug(b); printf("\n\n");
            abort();
        }

        arb_clear(a);
        arb_clear(b);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
