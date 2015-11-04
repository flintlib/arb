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

#include "arb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("contains_int....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        arb_t a;
        long c;
        int r, ok;

        arb_init(a);
        arb_randtest_special(a, state, 1 + n_randint(state, 500), 2);

        r = arb_contains_int(a);

        ok = !r;
        for (c = 0; c < 10; c++)
        {
            if (arb_contains_si(a, c) || arb_contains_si(a, -c))
            {
                ok = !ok;
                break;
            }
        }

        if (!ok)
        {
            printf("FAIL:\n\n");
            printf("a = "); arb_printd(a, 30); printf("\n\n");
            printf("r = %d\n\n", r);
            abort();
        }

        arb_clear(a);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
