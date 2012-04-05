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

    printf("set_fmpq....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        arb_t a;
        fmpq_t x;

        arb_init(a, n_randint(state, 200));
        fmpq_init(x);

        arb_randtest(a, state, 10);
        fmpq_randtest(x, state, 1 + n_randint(state, 200));
        arb_set_fmpq(a, x);

        if (!arb_contains_fmpq(a, x))
        {
            printf("FAIL\n\n");
            printf("a = "); arb_debug(a); printf("\n\n");
            printf("x = "); fmpq_print(x); printf("\n\n");
            abort();
        }

        arb_clear(a);
        fmpq_clear(x);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
