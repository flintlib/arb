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

#include "mag.h"
#include "long_extras.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("cmp....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        fmpr_t x;
        mag_t xb;
        long y;
        int c1, c2;

        fmpr_init(x);
        mag_init(xb);

        mag_randtest_special(xb, state, 100);
        y = z_randtest(state);

        mag_get_fmpr(x, xb);

        c1 = fmpr_cmp_2exp_si(x, y);
        c2 = mag_cmp_2exp_si(xb, y);

        if (c1 != c2)
        {
            printf("FAIL\n\n");
            printf("x = "); fmpr_print(x); printf("\n\n");
            printf("y = %ld", y);  printf("\n\n");
            printf("xb = "); mag_print(xb); printf("\n\n");
            abort();
        }

        fmpr_clear(x);
        mag_clear(xb);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

