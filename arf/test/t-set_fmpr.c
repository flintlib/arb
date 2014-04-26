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

    printf("set_fmpr....");
    fflush(stdout);

    flint_randinit(state);

    /* test exact roundtrip R -> Q -> R */
    for (iter = 0; iter < 100000; iter++)
    {
        long bits;
        arf_t x, z;
        fmpr_t y;

        bits = 2 + n_randint(state, 1000);

        arf_init(x);
        arf_init(z);
        fmpr_init(y);

        arf_randtest_special(x, state, bits, 1 + n_randint(state, 100));
        arf_randtest_special(z, state, bits, 1 + n_randint(state, 100));

        arf_get_fmpr(y, x);
        arf_set_fmpr(z, y);

        if (!arf_equal(x, z))
        {
            printf("FAIL\n\n");
            printf("bits: %ld\n", bits);
            printf("x = "); arf_print(x); printf("\n\n");
            printf("y = "); fmpr_print(y); printf("\n\n");
            printf("z = "); arf_print(z); printf("\n\n");
            abort();
        }

        arf_clear(x);
        arf_clear(z);
        fmpr_clear(y);
    }

    /* test exact roundtrip Q -> R -> Q */
    for (iter = 0; iter < 100000; iter++)
    {
        long bits;
        fmpr_t x, z;
        arf_t y;

        bits = 2 + n_randint(state, 1000);

        fmpr_init(x);
        fmpr_init(z);
        arf_init(y);

        fmpr_randtest_special(x, state, bits, 1 + n_randint(state, 100));
        fmpr_randtest_special(z, state, bits, 1 + n_randint(state, 100));

        arf_set_fmpr(y, x);
        arf_get_fmpr(z, y);

        if (!fmpr_equal(x, z))
        {
            printf("FAIL (2)\n\n");
            printf("bits: %ld\n", bits);
            printf("x = "); fmpr_print(x); printf("\n\n");
            printf("y = "); arf_print(y); printf("\n\n");
            printf("z = "); fmpr_print(z); printf("\n\n");
            abort();
        }

        fmpr_clear(x);
        fmpr_clear(z);
        arf_clear(y);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
