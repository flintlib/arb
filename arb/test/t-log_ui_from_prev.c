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

    printf("log_ui_from_prev....");
    fflush(stdout);
    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        arb_t z1, z2, z3;
        ulong n1, n2;
        long prec, accuracy;

        prec = 2 + n_randint(state, 3000);
        n1 = n_randint(state, 100000);
        n2 = n1 + 1 + n_randint(state, 1000);

        arb_init(z1);
        arb_init(z2);
        arb_init(z3);

        arb_log_ui(z1, n1, prec);
        arb_log_ui(z2, n2, prec);

        arb_log_ui_from_prev(z3, n2, z1, n1, prec);

        if (!arb_overlaps(z2, z3))
        {
            printf("FAIL: overlap\n\n");
            printf("prec = %ld, n1 = %lu, n2 = %lu\n\n", prec, n1, n2);
            printf("z1 = "); arb_printd(z1, prec / 3.33); printf("\n\n");
            printf("z2 = "); arb_printd(z2, prec / 3.33); printf("\n\n");
            printf("z3 = "); arb_printd(z3, prec / 3.33); printf("\n\n");
            abort();
        }

        accuracy = arb_rel_accuracy_bits(z3);

        if (accuracy < prec - 4)
        {
            printf("FAIL: accuracy = %ld, prec = %ld\n\n", accuracy, prec);
            printf("n1 = %lu, n2 = %lu\n\n", n1, n2);
            printf("z3 = "); arb_printd(z3, prec / 3.33); printf("\n\n");
            abort();
        }

        arb_clear(z1);
        arb_clear(z2);
        arb_clear(z3);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
