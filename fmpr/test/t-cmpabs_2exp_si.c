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

    printf("cmpabs_2exp_si....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        long bits, e;
        fmpr_t x, y;
        int cmp1, cmp2;

        bits = 2 + n_randint(state, 1000);
        e = n_randtest(state);

        fmpr_init(x);
        fmpr_init(y);

        fmpr_randtest_special(x, state, bits, 100);
        fmpr_set_ui_2exp_si(y, 1, e);

        cmp1 = fmpr_cmpabs(x, y);
        cmp2 = fmpr_cmpabs_2exp_si(x, e);

        if (cmp1 != cmp2)
        {
            printf("FAIL\n\n");
            printf("x = "); fmpr_print(x); printf("\n\n");
            printf("y = "); fmpr_print(y); printf("\n\n");
            printf("cmp1 = %d, cmp2 = %d\n\n", cmp1, cmp2);
            abort();
        }

        fmpr_clear(x);
        fmpr_clear(y);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

