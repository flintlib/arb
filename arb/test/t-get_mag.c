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

    printf("get_mag....");
    fflush(stdout);
    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        arb_t a, b;
        mag_t m;

        arb_init(a);
        arb_init(b);
        mag_init(m);

        arb_randtest_special(a, state, 200, 1 + n_randint(state, 100));
        arb_get_mag(m, a);
        MAG_CHECK_BITS(m)

        if (arf_is_nan(arb_midref(a)))
            arf_nan(arb_midref(b));
        else
            arf_zero(arb_midref(b));
        mag_set(arb_radref(b), m);

        if (!arb_contains(b, a))
        {
            printf("FAIL:\n\n");
            printf("a = "); arb_print(a); printf("\n\n");
            printf("b = "); arb_print(b); printf("\n\n");
            printf("m = "); mag_print(m); printf("\n\n");
            abort();
        }

        arb_clear(a);
        arb_clear(b);
        mag_clear(m);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

