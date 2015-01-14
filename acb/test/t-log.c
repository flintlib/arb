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

    printf("log....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        acb_t x, a, b;
        long prec1, prec2;

        prec1 = 2 + n_randint(state, 1000);
        prec2 = prec1 + 30;

        acb_init(x);
        acb_init(a);
        acb_init(b);

        acb_randtest_special(x, state, 1 + n_randint(state, 1000), 2 + n_randint(state, 100));
        acb_randtest_special(a, state, 1 + n_randint(state, 1000), 2 + n_randint(state, 100));
        acb_randtest_special(b, state, 1 + n_randint(state, 1000), 2 + n_randint(state, 100));

        acb_log(a, x, prec1);
        acb_log(b, x, prec2);

        /* check consistency */
        if (!acb_overlaps(a, b))
        {
            printf("FAIL: overlap\n\n");
            printf("x = "); acb_printd(x, 15); printf("\n\n");
            printf("a = "); acb_printd(a, 15); printf("\n\n");
            printf("b = "); acb_printd(b, 15); printf("\n\n");
            abort();
        }

        /* check exp(log(x)) = x */
        acb_exp(b, b, prec1);

        if (!acb_contains(b, x))
        {
            printf("FAIL: functional equation\n\n");
            printf("x = "); acb_printd(x, 15); printf("\n\n");
            printf("b = "); acb_printd(b, 15); printf("\n\n");
            abort();
        }

        acb_log(x, x, prec1);

        if (!acb_overlaps(a, x))
        {
            printf("FAIL: aliasing\n\n");
            printf("a = "); acb_printd(a, 15); printf("\n\n");
            printf("x = "); acb_printd(x, 15); printf("\n\n");
            abort();
        }

        acb_clear(x);
        acb_clear(a);
        acb_clear(b);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

