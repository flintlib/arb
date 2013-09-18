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

#include "fmprb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("trim....");
    fflush(stdout);
    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        fmprb_t x, y;
        long acc1, acc2;
        int accuracy_ok;

        fmprb_init(x);
        fmprb_init(y);

        fmprb_randtest_special(x, state, 1000, 100);

        fmprb_trim(y, x);

        if (!fmprb_contains(y, x))
        {
            printf("FAIL (containment):\n\n");
            printf("x = "); fmprb_print(x); printf("\n\n");
            printf("y = "); fmprb_print(y); printf("\n\n");
            abort();
        }

        acc1 = fmprb_rel_accuracy_bits(x);
        acc2 = fmprb_rel_accuracy_bits(y);

        accuracy_ok = (acc1 < 0 && acc2 <= acc1) || (acc2 >= acc1 - 1);

        if (!accuracy_ok)
        {
            printf("FAIL (accuracy):\n\n");
            printf("x: %ld, y = %ld\n\n", acc1, acc2);
            printf("x = "); fmprb_print(x); printf("\n\n");
            printf("y = "); fmprb_print(y); printf("\n\n");
            abort();
        }

        fmprb_trim(x, x);

        if (!fmprb_equal(y, x))
        {
            printf("FAIL (aliasing):\n\n");
            printf("x = "); fmprb_print(x); printf("\n\n");
            printf("y = "); fmprb_print(y); printf("\n\n");
            abort();
        }

        fmprb_clear(x);
        fmprb_clear(y);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

