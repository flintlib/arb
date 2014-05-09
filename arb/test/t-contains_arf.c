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

    printf("contains_arf....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        arb_t a;
        arf_t b;
        fmpq_t am, ar, bm, t;
        int c1, c2;

        arb_init(a);
        arf_init(b);

        fmpq_init(am);
        fmpq_init(ar);
        fmpq_init(bm);
        fmpq_init(t);

        arb_randtest(a, state, 1 + n_randint(state, 500), 14);
        arf_randtest(b, state, 1 + n_randint(state, 500), 14);

        arf_get_fmpq(am, arb_midref(a));
        mag_get_fmpq(ar, arb_radref(a));
        arf_get_fmpq(bm, b);

        c1 = arb_contains_arf(a, b);

        fmpq_sub(t, am, ar);
        c2 = fmpq_cmp(t, bm) <= 0;

        fmpq_add(t, am, ar);
        c2 = c2 && (fmpq_cmp(t, bm) >= 0);

        if (c1 != c2)
        {
            printf("FAIL:\n\n");
            printf("a = "); arb_print(a); printf("\n\n");
            printf("b = "); arf_print(b); printf("\n\n");
            printf("am = "); fmpq_print(am); printf("\n\n");
            printf("ar = "); fmpq_print(ar); printf("\n\n");
            printf("bm = "); fmpq_print(bm); printf("\n\n");
            printf("t = "); fmpq_print(t); printf("\n\n");
            printf("c1 = %d, c2 = %d\n\n", c1, c2);
            abort();
        }

        arb_clear(a);
        arf_clear(b);

        fmpq_clear(am);
        fmpq_clear(ar);
        fmpq_clear(bm);
        fmpq_clear(t);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
