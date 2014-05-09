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

#include "arb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("get_rand_fmpq....");
    fflush(stdout);
    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        arb_t x;
        fmpq_t q;

        arb_init(x);
        fmpq_init(q);

        arb_randtest(x, state, 1 + n_randint(state, 200), 10);
        arb_get_rand_fmpq(q, state, x, 1 + n_randint(state, 200));

        if (!arb_contains_fmpq(x, q) || !fmpq_is_canonical(q))
        {
            printf("FAIL:\n\n");
            printf("x = "); arb_print(x); printf("\n\n");
            printf("q = "); fmpq_print(q); printf("\n\n");
            abort();
        }

        arb_clear(x);
        fmpq_clear(q);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
