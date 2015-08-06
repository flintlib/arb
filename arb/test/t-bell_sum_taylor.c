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

    Copyright (C) 2015 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include "arb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("bell_sum_taylor....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        arb_t s1, s2;
        fmpz_t a, b, n;
        long prec;

        arb_init(s1);
        arb_init(s2);
        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(n);

        prec = 2 + n_randint(state, 300);
        fmpz_randtest_unsigned(n, state, 1 + n_randint(state, 100));
        fmpz_randtest_unsigned(a, state, 1 + n_randint(state, 100));
        fmpz_add_ui(b, a, n_randint(state, 100));

        arb_bell_sum_bsplit(s1, n, a, b, NULL, prec);
        arb_bell_sum_taylor(s2, n, a, b, NULL, prec);

        if (!arb_overlaps(s1, s2) || (arb_rel_accuracy_bits(s1) < prec - 4)
            || (arb_rel_accuracy_bits(s2) < prec - 4))
        {
            printf("FAIL: overlap or accuracy\n\n");
            printf("prec = %ld\n\n", prec);
            printf("n = "); fmpz_print(n); printf("\n\n");
            printf("a = "); fmpz_print(a); printf("\n\n");
            printf("b = "); fmpz_print(b); printf("\n\n");
            printf("s1 = "); arb_printn(s1, 100, 0); printf("\n\n");
            printf("s2 = "); arb_printn(s2, 100, 0); printf("\n\n");
            abort();
        }

        arb_clear(s1);
        arb_clear(s2);
        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(n);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

