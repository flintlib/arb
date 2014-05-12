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

#include "arb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("coth....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        arb_t x, y, a, b, c, d;
        long prec1, prec2;

        prec1 = 2 + n_randint(state, 1000);
        prec2 = prec1 + 30;

        arb_init(x);
        arb_init(y);
        arb_init(a);
        arb_init(b);
        arb_init(c);
        arb_init(d);

        arb_randtest_precise(x, state, 1 + n_randint(state, 1000), 100);
        arb_randtest_precise(y, state, 1 + n_randint(state, 1000), 100);

        arb_coth(a, x, prec1);
        arb_coth(b, x, prec2);

        /* check consistency */
        if (!arb_overlaps(a, b))
        {
            printf("FAIL: overlap\n\n");
            printf("x = "); arb_print(x); printf("\n\n");
            printf("a = "); arb_print(a); printf("\n\n");
            printf("b = "); arb_print(b); printf("\n\n");
            abort();
        }

        /* check coth(x+y) = (1 + coth(x) coth(y)) / (coth(x) + coth(y)) */
        arb_add(b, x, y, prec1);
        arb_coth(b, b, prec1);

        arb_coth(c, y, prec1);
        arb_add(d, a, c, prec1);
        arb_mul(c, a, c, prec1);
        arb_add_ui(c, c, 1, prec1);
        arb_div(d, c, d, prec1);

        if (!arb_overlaps(b, d))
        {
            printf("FAIL: functional equation\n\n");
            printf("x = "); arb_print(x); printf("\n\n");
            printf("y = "); arb_print(y); printf("\n\n");
            printf("b = "); arb_print(b); printf("\n\n");
            printf("d = "); arb_print(d); printf("\n\n");
            abort();
        }

        arb_coth(x, x, prec1);

        if (!arb_overlaps(a, x))
        {
            printf("FAIL: aliasing\n\n");
            printf("a = "); arb_print(a); printf("\n\n");
            printf("x = "); arb_print(x); printf("\n\n");
            abort();
        }

        arb_clear(x);
        arb_clear(y);
        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        arb_clear(d);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

