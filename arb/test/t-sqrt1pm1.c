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

    Copyright (C) 2012, 2013 Fredrik Johansson

******************************************************************************/

#include "arb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("sqrt1pm1....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 20000; iter++)
    {
        arb_t a, b, c, d;
        long prec0, prec1, prec2;

        if (iter % 10 == 0)
            prec0 = 10000;
        else
            prec0 = 1000;

        prec1 = 2 + n_randint(state, prec0);
        prec2 = 2 + n_randint(state, prec0);

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arb_init(d);

        arb_randtest_special(a, state, 1 + n_randint(state, prec0), 100);
        arb_randtest_special(b, state, 1 + n_randint(state, prec0), 100);
        arb_randtest_special(c, state, 1 + n_randint(state, prec0), 100);

        arb_sqrt1pm1(b, a, prec1);
        arb_sqrt1pm1(c, a, prec2);

        if (!arb_overlaps(b, c))
        {
            printf("FAIL: overlap\n\n");
            printf("a = "); arb_print(a); printf("\n\n");
            printf("b = "); arb_print(b); printf("\n\n");
            printf("c = "); arb_print(c); printf("\n\n");
            abort();
        }

        /* compare with sqrt */
        arb_add_ui(d, a, 1, prec2);
        arb_sqrt(d, d, prec2);
        arb_sub_ui(d, d, 1, prec2);

        if (!arb_overlaps(c, d))
        {
            printf("FAIL: comparison with log\n\n");
            printf("a = "); arb_print(a); printf("\n\n");
            printf("b = "); arb_print(b); printf("\n\n");
            printf("c = "); arb_print(c); printf("\n\n");
            printf("d = "); arb_print(d); printf("\n\n");
            abort();
        }

        arb_sqrt1pm1(a, a, prec1);

        if (!arb_overlaps(a, b))
        {
            printf("FAIL: aliasing\n\n");
            printf("a = "); arb_print(a); printf("\n\n");
            printf("b = "); arb_print(b); printf("\n\n");
            abort();
        }

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
