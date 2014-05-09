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

    printf("union....");
    fflush(stdout);
    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        arb_t x, y, z;
        long prec;
        int alias;

        arb_init(x);
        arb_init(y);
        arb_init(z);

        arb_randtest_special(x, state, 200, 10);
        arb_randtest_special(y, state, 200, 10);
        arb_randtest_special(z, state, 200, 10);

        prec = 2 + n_randint(state, 200);

        arb_union(z, x, y, prec);

        if (!arb_contains(z, x) || !arb_contains(z, y))
        {
            printf("FAIL:\n\n");
            printf("x = "); arb_print(x); printf("\n\n");
            printf("y = "); arb_print(y); printf("\n\n");
            printf("z = "); arb_print(z); printf("\n\n");
            abort();
        }

        if (n_randint(state, 2))
        {
            arb_union(x, x, y, prec);
            alias = arb_equal(x, z);
        }
        else
        {
            arb_union(y, x, y, prec);
            alias = arb_equal(y, z);
        }

        if (!alias)
        {
            printf("FAIL (aliasing):\n\n");
            printf("x = "); arb_print(x); printf("\n\n");
            printf("y = "); arb_print(y); printf("\n\n");
            printf("z = "); arb_print(z); printf("\n\n");
            abort();
        }

        arb_clear(x);
        arb_clear(y);
        arb_clear(z);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

