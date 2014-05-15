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

#include "acb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("tan....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        acb_t x, y, a, b, c, d;
        long prec1, prec2;

        prec1 = 2 + n_randint(state, 1000);
        prec2 = prec1 + 30;

        acb_init(x);
        acb_init(y);
        acb_init(a);
        acb_init(b);
        acb_init(c);
        acb_init(d);

        acb_randtest(x, state, 1 + n_randint(state, 1000), 2 + n_randint(state, 100));
        acb_randtest(y, state, 1 + n_randint(state, 1000), 2 + n_randint(state, 100));

        acb_tan(a, x, prec1);
        acb_tan(b, x, prec2);

        /* check consistency */
        if (!acb_overlaps(a, b))
        {
            printf("FAIL: overlap\n\n");
            printf("x = "); acb_print(x); printf("\n\n");
            printf("a = "); acb_print(a); printf("\n\n");
            printf("b = "); acb_print(b); printf("\n\n");
            abort();
        }

        /* check tan(x+y) = (tan(x) + tan(y)) / (1 - tan(x) tan(y)) */
        acb_add(b, x, y, prec1);
        acb_tan(b, b, prec1);

        acb_tan(c, y, prec1);
        acb_add(d, a, c, prec1);
        acb_mul(c, a, c, prec1);
        acb_sub_ui(c, c, 1, prec1);
        acb_neg(c, c);
        acb_div(d, d, c, prec1);

        if (!acb_overlaps(b, d))
        {
            printf("FAIL: functional equation\n\n");
            printf("x = "); acb_print(x); printf("\n\n");
            printf("y = "); acb_print(y); printf("\n\n");
            printf("b = "); acb_print(b); printf("\n\n");
            printf("d = "); acb_print(d); printf("\n\n");
            abort();
        }

        acb_tan(x, x, prec1);

        if (!acb_overlaps(a, x))
        {
            printf("FAIL: aliasing\n\n");
            printf("a = "); acb_print(a); printf("\n\n");
            printf("x = "); acb_print(x); printf("\n\n");
            abort();
        }

        acb_clear(x);
        acb_clear(y);
        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
        acb_clear(d);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

