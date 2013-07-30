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

#include "fmprb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("coth....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        fmprb_t x, y, a, b, c, d;
        long prec1, prec2;

        prec1 = 2 + n_randint(state, 1000);
        prec2 = prec1 + 30;

        fmprb_init(x);
        fmprb_init(y);
        fmprb_init(a);
        fmprb_init(b);
        fmprb_init(c);
        fmprb_init(d);

        fmprb_randtest_precise(x, state, 1 + n_randint(state, 1000), 100);
        fmprb_randtest_precise(y, state, 1 + n_randint(state, 1000), 100);

        fmprb_coth(a, x, prec1);
        fmprb_coth(b, x, prec2);

        /* check consistency */
        if (!fmprb_overlaps(a, b))
        {
            printf("FAIL: overlap\n\n");
            printf("x = "); fmprb_print(x); printf("\n\n");
            printf("a = "); fmprb_print(a); printf("\n\n");
            printf("b = "); fmprb_print(b); printf("\n\n");
            abort();
        }

        /* check coth(x+y) = (1 + coth(x) coth(y)) / (coth(x) + coth(y)) */
        fmprb_add(b, x, y, prec1);
        fmprb_coth(b, b, prec1);

        fmprb_coth(c, y, prec1);
        fmprb_add(d, a, c, prec1);
        fmprb_mul(c, a, c, prec1);
        fmprb_add_ui(c, c, 1, prec1);
        fmprb_div(d, c, d, prec1);

        if (!fmprb_overlaps(b, d))
        {
            printf("FAIL: functional equation\n\n");
            printf("x = "); fmprb_print(x); printf("\n\n");
            printf("y = "); fmprb_print(y); printf("\n\n");
            printf("b = "); fmprb_print(b); printf("\n\n");
            printf("d = "); fmprb_print(d); printf("\n\n");
            abort();
        }

        fmprb_coth(x, x, prec1);

        if (!fmprb_overlaps(a, x))
        {
            printf("FAIL: aliasing\n\n");
            printf("a = "); fmprb_print(a); printf("\n\n");
            printf("x = "); fmprb_print(x); printf("\n\n");
            abort();
        }

        fmprb_clear(x);
        fmprb_clear(y);
        fmprb_clear(a);
        fmprb_clear(b);
        fmprb_clear(c);
        fmprb_clear(d);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

