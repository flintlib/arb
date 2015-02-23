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

#include "acb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("rel_accuracy_bits....");
    fflush(stdout);

    flint_randinit(state);

    /* test aliasing of c and a */
    for (iter = 0; iter < 10000; iter++)
    {
        arb_t x;
        acb_t z;
        long a1, a2;

        arb_init(x);
        acb_init(z);

        arb_randtest_special(x, state, 1 + n_randint(state, 200), 1 + n_randint(state, 200));
        acb_set_arb(z, x);

        a1 = arb_rel_accuracy_bits(x);
        a2 = acb_rel_accuracy_bits(z);

        if (a1 != a2)
        {
            printf("FAIL: acb != arb\n\n");
            printf("x = "); arb_print(x); printf("\n\n");
            printf("z = "); acb_print(z); printf("\n\n");
            printf("a1 = %ld, a2 = %ld\n\n", a1, a2);
            abort();
        }

        acb_randtest_special(z, state, 1 + n_randint(state, 200), 1 + n_randint(state, 200));

        a1 = acb_rel_accuracy_bits(z);

        if (n_randint(state, 2))
            arf_swap(arb_midref(acb_realref(z)), arb_midref(acb_imagref(z)));

        if (n_randint(state, 2))
            mag_swap(arb_radref(acb_realref(z)), arb_radref(acb_imagref(z)));

        a2 = acb_rel_accuracy_bits(z);

        if (a1 != a2)
        {
            printf("FAIL: swapping\n\n");
            printf("z = "); acb_print(z); printf("\n\n");
            printf("a1 = %ld, a2 = %ld\n\n", a1, a2);
            abort();
        }

        acb_randtest_special(z, state, 1 + n_randint(state, 200), 1 + n_randint(state, 200));

        if (arf_cmpabs(arb_midref(acb_realref(z)), arb_midref(acb_imagref(z))) >= 0)
            arf_set(arb_midref(x), arb_midref(acb_realref(z)));
        else
            arf_set(arb_midref(x), arb_midref(acb_imagref(z)));

        if (mag_cmp(arb_radref(acb_realref(z)), arb_radref(acb_imagref(z))) >= 0)
            mag_set(arb_radref(x), arb_radref(acb_realref(z)));
        else
            mag_set(arb_radref(x), arb_radref(acb_imagref(z)));

        a1 = acb_rel_accuracy_bits(z);
        a2 = arb_rel_accuracy_bits(x);

        if (a1 != a2)
        {
            printf("FAIL: acb != arb (2)\n\n");
            printf("x = "); arb_print(x); printf("\n\n");
            printf("z = "); acb_print(z); printf("\n\n");
            printf("a1 = %ld, a2 = %ld\n\n", a1, a2);
            abort();
        }

        arb_clear(x);
        acb_clear(z);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

