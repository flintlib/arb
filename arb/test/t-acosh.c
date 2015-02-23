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

    printf("acosh....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        arb_t x, a, b;
        long prec1, prec2;

        prec1 = 2 + n_randint(state, 1000);
        prec2 = prec1 + 30;

        arb_init(x);
        arb_init(a);
        arb_init(b);

        arb_randtest_special(x, state, 1 + n_randint(state, 1000), 2 + n_randint(state, 100));
        arb_randtest_special(a, state, 1 + n_randint(state, 1000), 2 + n_randint(state, 100));
        arb_randtest_special(b, state, 1 + n_randint(state, 1000), 2 + n_randint(state, 100));

        arb_acosh(a, x, prec1);
        arb_acosh(b, x, prec2);

        /* check consistency */
        if (!arb_overlaps(a, b))
        {
            printf("FAIL: overlap\n\n");
            printf("x = "); arb_printd(x, 15); printf("\n\n");
            printf("a = "); arb_printd(a, 15); printf("\n\n");
            printf("b = "); arb_printd(b, 15); printf("\n\n");
            abort();
        }

        /* check cosh(acosh(x)) = x */
        arb_cosh(b, b, prec1);

        if (!arb_contains(b, x))
        {
            printf("FAIL: functional equation\n\n");
            printf("x = "); arb_printd(x, 15); printf("\n\n");
            printf("b = "); arb_printd(b, 15); printf("\n\n");
            abort();
        }

        arb_acosh(x, x, prec1);

        if (!arb_overlaps(a, x))
        {
            printf("FAIL: aliasing\n\n");
            printf("a = "); arb_printd(a, 15); printf("\n\n");
            printf("x = "); arb_printd(x, 15); printf("\n\n");
            abort();
        }

        arb_clear(x);
        arb_clear(a);
        arb_clear(b);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

