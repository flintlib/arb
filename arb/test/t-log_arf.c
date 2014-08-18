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

#include "arb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("log_arf....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 5000; iter++)
    {
        arf_t x;
        arb_t y1, y2;
        long prec1, prec2, acc1, acc2;

        prec1 = 2 + n_randint(state, 9000);
        prec2 = 2 + n_randint(state, 9000);

        arf_init(x);
        arb_init(y1);
        arb_init(y2);

        arf_randtest_special(x, state, 1 + n_randint(state, 9000), 200);
        arb_randtest_special(y1, state, 1 + n_randint(state, 9000), 200);
        arb_randtest_special(y2, state, 1 + n_randint(state, 9000), 200);

        if (n_randint(state, 2))
            arf_add_ui(x, x, 1, 2 + n_randint(state, 9000), ARF_RND_DOWN);

        arb_log_arf(y1, x, prec1);
        arb_log_arf(y2, x, prec2);

        if (!arb_overlaps(y1, y2))
        {
            printf("FAIL: overlap\n\n");
            printf("prec1 = %ld, prec2 = %ld\n\n", prec1, prec2);
            printf("x = "); arf_print(x); printf("\n\n");
            printf("y1 = "); arb_print(y1); printf("\n\n");
            printf("y2 = "); arb_print(y2); printf("\n\n");
            abort();
        }

        acc1 = arb_rel_accuracy_bits(y1);
        acc2 = arb_rel_accuracy_bits(y2);

        if (arf_sgn(x) > 0)
        {
            if (acc1 < prec1 - 2 || acc2 < prec2 - 2)
            {
                printf("FAIL: accuracy\n\n");
                printf("prec1 = %ld, prec2 = %ld\n\n", prec1, prec2);
                printf("acc1 = %ld, acc2 = %ld\n\n", acc1, acc2);
                printf("x = "); arf_print(x); printf("\n\n");
                printf("y1 = "); arb_print(y1); printf("\n\n");
                printf("y2 = "); arb_print(y2); printf("\n\n");
                abort();
            }
        }

        arf_clear(x);
        arb_clear(y1);
        arb_clear(y2);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

