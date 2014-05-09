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

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include "arb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("trim....");
    fflush(stdout);
    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        arb_t x, y;
        long acc1, acc2;
        int accuracy_ok;

        arb_init(x);
        arb_init(y);

        arb_randtest_special(x, state, 1000, 100);

        arb_trim(y, x);

        if (!arb_contains(y, x))
        {
            printf("FAIL (containment):\n\n");
            printf("x = "); arb_print(x); printf("\n\n");
            printf("y = "); arb_print(y); printf("\n\n");
            abort();
        }

        acc1 = arb_rel_accuracy_bits(x);
        acc2 = arb_rel_accuracy_bits(y);

        accuracy_ok = (acc1 < 0 && acc2 <= acc1) || (acc2 >= acc1 - 1);

        if (!accuracy_ok)
        {
            printf("FAIL (accuracy):\n\n");
            printf("x: %ld, y = %ld\n\n", acc1, acc2);
            printf("x = "); arb_print(x); printf("\n\n");
            printf("y = "); arb_print(y); printf("\n\n");
            abort();
        }

        arb_trim(x, x);

        if (!arb_equal(y, x))
        {
            printf("FAIL (aliasing):\n\n");
            printf("x = "); arb_print(x); printf("\n\n");
            printf("y = "); arb_print(y); printf("\n\n");
            abort();
        }

        arb_clear(x);
        arb_clear(y);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

