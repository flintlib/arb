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

#include <stdio.h>
#include "arb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("bernoulli_ui....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        arb_t b1, b2;
        ulong n;
        long prec1, prec2, acc1, acc2;

        n = n_randint(state, 10000);
        prec1 = 2 + n_randint(state, 10000);
        prec2 = prec1 + 100;

        arb_init(b1);
        arb_init(b2);

        arb_bernoulli_ui(b1, n, prec1);
        arb_bernoulli_ui(b2, n, prec2);

        if (!arb_overlaps(b1, b2))
        {
            printf("FAIL: overlap\n\n");
            printf("n = %lu\n\n", n);
            printf("b1 = "); arb_print(b1); printf("\n\n");
            printf("b2 = "); arb_print(b2); printf("\n\n");
            abort();
        }

        acc1 = arb_rel_accuracy_bits(b1);
        acc2 = arb_rel_accuracy_bits(b2);

        if (acc1 < prec1 - 2 || acc2 < prec2 - 2)
        {
            printf("FAIL: poor accuracy\n\n");
            printf("prec1 = %ld\n", prec1);
            printf("prec2 = %ld\n", prec2);
            printf("b1 = "); arb_printd(b1, prec1 / 3.33); printf("\n\n");
            printf("b2 = "); arb_printd(b2, prec2 / 3.33); printf("\n\n");
            abort();
        }

        arb_clear(b1);
        arb_clear(b2);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

