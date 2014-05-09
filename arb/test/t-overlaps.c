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

    printf("overlaps....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        arb_t a, b;
        fmpq_t am, ar, bm, br, t, u;
        int c1, c2;

        arb_init(a);
        arb_init(b);

        fmpq_init(am);
        fmpq_init(ar);
        fmpq_init(bm);
        fmpq_init(br);
        fmpq_init(t);
        fmpq_init(u);

        arb_randtest(a, state, 1 + n_randint(state, 500), 14);
        arb_randtest(b, state, 1 + n_randint(state, 500), 14);

        arf_get_fmpq(am, arb_midref(a));
        mag_get_fmpq(ar, arb_radref(a));
        arf_get_fmpq(bm, arb_midref(b));
        mag_get_fmpq(br, arb_radref(b));

        fmpq_sub(t, am, bm);
        fmpz_abs(fmpq_numref(t), fmpq_numref(t));
        fmpq_add(u, ar, br);

        c1 = arb_overlaps(a, b);

        c2 = (fmpq_cmp(t, u) <= 0);

        if (c1 != c2)
        {
            printf("FAIL:\n\n");
            printf("a = "); arb_print(a); printf("\n\n");
            printf("b = "); arb_print(b); printf("\n\n");
            printf("am = "); fmpq_print(am); printf("\n\n");
            printf("ar = "); fmpq_print(ar); printf("\n\n");
            printf("bm = "); fmpq_print(bm); printf("\n\n");
            printf("br = "); fmpq_print(br); printf("\n\n");
            printf("t = "); fmpq_print(t); printf("\n\n");
            printf("u = "); fmpq_print(u); printf("\n\n");
            printf("c1 = %d, c2 = %d\n\n", c1, c2);
            abort();
        }

        arb_clear(a);
        arb_clear(b);

        fmpq_clear(am);
        fmpq_clear(ar);
        fmpq_clear(bm);
        fmpq_clear(br);
        fmpq_clear(t);
        fmpq_clear(u);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
