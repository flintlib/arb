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
    slong iter;
    flint_rand_t state;

    flint_printf("trim....");
    fflush(stdout);
    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        arb_t x, y;
        slong acc1, acc2;
        int accuracy_ok;

        arb_init(x);
        arb_init(y);

        arb_randtest_special(x, state, 1000, 100);

        arb_trim(y, x);

        if (!arb_contains(y, x))
        {
            flint_printf("FAIL (containment):\n\n");
            flint_printf("x = "); arb_print(x); flint_printf("\n\n");
            flint_printf("y = "); arb_print(y); flint_printf("\n\n");
            abort();
        }

        acc1 = arb_rel_accuracy_bits(x);
        acc2 = arb_rel_accuracy_bits(y);

        accuracy_ok = (acc1 < 0 && acc2 <= acc1) || (acc2 >= acc1 - 1);

        if (!accuracy_ok)
        {
            flint_printf("FAIL (accuracy):\n\n");
            flint_printf("x: %wd, y = %wd\n\n", acc1, acc2);
            flint_printf("x = "); arb_print(x); flint_printf("\n\n");
            flint_printf("y = "); arb_print(y); flint_printf("\n\n");
            abort();
        }

        arb_trim(x, x);

        if (!arb_equal(y, x))
        {
            flint_printf("FAIL (aliasing):\n\n");
            flint_printf("x = "); arb_print(x); flint_printf("\n\n");
            flint_printf("y = "); arb_print(y); flint_printf("\n\n");
            abort();
        }

        arb_clear(x);
        arb_clear(y);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

