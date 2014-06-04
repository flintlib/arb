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

#include "arf.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("get_d....");
    fflush(stdout);

    flint_randinit(state);

    /* test exact roundtrip */
    for (iter = 0; iter < 100000; iter++)
    {
        arf_t x, z;
        double y;

        arf_init(x);
        arf_init(z);

        arf_randtest_special(x, state, 53, 8);
        y = arf_get_d(x, ARF_RND_DOWN);
        arf_set_d(z, y);

        if (!arf_equal(x, z))
        {
            printf("FAIL:\n\n");
            printf("x = "); arf_print(x); printf("\n\n");
            printf("y = %.17g\n\n", y);
            printf("z = "); arf_print(z); printf("\n\n");
            abort();
        }

        arf_clear(x);
        arf_clear(z);
    }

    /* test rounding */
    for (iter = 0; iter < 100000; iter++)
    {
        arf_t x, z, w;
        arf_rnd_t rnd;
        double y;

        arf_init(x);
        arf_init(z);
        arf_init(w);

        arf_randtest_special(x, state, 300, 8);

        switch (n_randint(state, 4))
        {
            case 0:  rnd = ARF_RND_DOWN; break;
            case 1:  rnd = ARF_RND_UP; break;
            case 2:  rnd = ARF_RND_FLOOR; break;
            default: rnd = ARF_RND_CEIL; break;
        }

        y = arf_get_d(x, rnd);
        arf_set_d(w, y);

        arf_set_round(z, x, 53, rnd);

        if (!arf_equal(w, z))
        {
            printf("FAIL:\n\n");
            printf("x = "); arf_print(x); printf("\n\n");
            printf("y = %.17g\n\n", y);
            printf("z = "); arf_print(z); printf("\n\n");
            printf("w = "); arf_print(w); printf("\n\n");
            abort();
        }

        arf_clear(x);
        arf_clear(z);
        arf_clear(w);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

