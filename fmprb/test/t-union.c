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
    slong iter;
    flint_rand_t state;

    flint_printf("union....");
    fflush(stdout);
    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        fmprb_t x, y, z;
        slong prec;
        int alias;

        fmprb_init(x);
        fmprb_init(y);
        fmprb_init(z);

        fmprb_randtest_special(x, state, 200, 10);
        fmprb_randtest_special(y, state, 200, 10);
        fmprb_randtest_special(z, state, 200, 10);

        prec = 2 + n_randint(state, 200);

        fmprb_union(z, x, y, prec);

        if (!fmprb_contains(z, x) || !fmprb_contains(z, y))
        {
            flint_printf("FAIL:\n\n");
            flint_printf("x = "); fmprb_print(x); flint_printf("\n\n");
            flint_printf("y = "); fmprb_print(y); flint_printf("\n\n");
            flint_printf("z = "); fmprb_print(z); flint_printf("\n\n");
            abort();
        }

        if (n_randint(state, 2))
        {
            fmprb_union(x, x, y, prec);
            alias = fmprb_equal(x, z);
        }
        else
        {
            fmprb_union(y, x, y, prec);
            alias = fmprb_equal(y, z);
        }

        if (!alias)
        {
            flint_printf("FAIL (aliasing):\n\n");
            flint_printf("x = "); fmprb_print(x); flint_printf("\n\n");
            flint_printf("y = "); fmprb_print(y); flint_printf("\n\n");
            flint_printf("z = "); fmprb_print(z); flint_printf("\n\n");
            abort();
        }

        fmprb_clear(x);
        fmprb_clear(y);
        fmprb_clear(z);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

