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

#include "fmprb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("set_interval_fmpr....");
    fflush(stdout);
    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        fmprb_t x;
        fmpr_t a, b;

        fmprb_init(x);
        fmpr_init(a);
        fmpr_init(b);

        fmpr_randtest_special(a, state, 200, 10);
        fmpr_randtest_special(b, state, 200, 10);
        if (fmpr_cmp(a, b) > 0)
            fmpr_swap(a, b);

        fmprb_set_interval_fmpr(x, a, b, 2 + n_randint(state, 200));

        if (!fmprb_contains_fmpr(x, a) || !fmprb_contains_fmpr(x, b))
        {
            printf("FAIL:\n\n");
            printf("x = "); fmprb_print(x); printf("\n\n");
            printf("a = "); fmpr_print(a); printf("\n\n");
            printf("b = "); fmpr_print(b); printf("\n\n");
            abort();
        }

        fmprb_clear(x);
        fmpr_clear(a);
        fmpr_clear(b);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

