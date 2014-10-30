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

    Copyright (C) 2014 Fredrik Johansson

******************************************************************************/

#include "arb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("ceil....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        arb_t a, b;
        fmpq_t x;
        fmpz_t y;
        long prec;

        arb_init(a);
        arb_init(b);

        fmpq_init(x);
        fmpz_init(y);

        arb_randtest(a, state, 1 + n_randint(state, 200), 10);
        arb_randtest(b, state, 1 + n_randint(state, 200), 10);
        prec = 2 + n_randint(state, 200);

        arb_get_rand_fmpq(x, state, a, 1 + n_randint(state, 200));

        arb_ceil(b, a, prec);
        fmpz_cdiv_q(y, fmpq_numref(x), fmpq_denref(x));

        if (!arb_contains_fmpz(b, y))
        {
            printf("FAIL: containment\n\n");
            printf("a = "); arb_print(a); printf("\n\n");
            printf("x = "); fmpq_print(x); printf("\n\n");
            printf("b = "); arb_print(b); printf("\n\n");
            printf("y = "); fmpz_print(y); printf("\n\n");
            abort();
        }

        arb_ceil(a, a, prec);

        if (!arb_equal(a, b))
        {
            printf("FAIL: aliasing\n\n");
            printf("a = "); arb_print(a); printf("\n\n");
            printf("b = "); arb_print(b); printf("\n\n");
            abort();
        }

        arb_clear(a);
        arb_clear(b);

        fmpq_clear(x);
        fmpz_clear(y);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
