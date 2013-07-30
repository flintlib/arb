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

#include "fmpr.h"

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
        fmpr_t x, z;
        double y;

        fmpr_init(x);
        fmpr_init(z);

        fmpr_randtest_special(x, state, 53, 8);
        y = fmpr_get_d(x, FMPR_RND_DOWN);
        fmpr_set_d(z, y);

        if (!fmpr_equal(x, z))
        {
            printf("FAIL:\n\n");
            printf("x = "); fmpr_print(x); printf("\n\n");
            printf("y = %.17g\n\n", y);
            printf("z = "); fmpr_print(z); printf("\n\n");
            abort();
        }

        fmpr_clear(x);
        fmpr_clear(z);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

