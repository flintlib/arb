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

#include "fmprb.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("contains_fmpq....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        fmprb_t a;
        fmpr_t b;
        fmpq_t am, ar, bm, t;
        int c1, c2;

        fmprb_init(a);
        fmpr_init(b);

        fmpq_init(am);
        fmpq_init(ar);
        fmpq_init(bm);
        fmpq_init(t);

        fmprb_randtest(a, state, 1 + n_randint(state, 500), 14);
        fmpr_randtest(b, state, 1 + n_randint(state, 500), 14);

        fmpr_get_fmpq(am, fmprb_midref(a));
        fmpr_get_fmpq(ar, fmprb_radref(a));
        fmpr_get_fmpq(bm, b);

        c1 = fmprb_contains_fmpr(a, b);

        fmpq_sub(t, am, ar);
        c2 = fmpq_cmp(t, bm) <= 0;

        fmpq_add(t, am, ar);
        c2 = c2 && (fmpq_cmp(t, bm) >= 0);

        if (c1 != c2)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "); fmprb_print(a); flint_printf("\n\n");
            flint_printf("b = "); fmpr_print(b); flint_printf("\n\n");
            flint_printf("am = "); fmpq_print(am); flint_printf("\n\n");
            flint_printf("ar = "); fmpq_print(ar); flint_printf("\n\n");
            flint_printf("bm = "); fmpq_print(bm); flint_printf("\n\n");
            flint_printf("t = "); fmpq_print(t); flint_printf("\n\n");
            flint_printf("c1 = %d, c2 = %d\n\n", c1, c2);
            abort();
        }

        fmprb_clear(a);
        fmpr_clear(b);

        fmpq_clear(am);
        fmpq_clear(ar);
        fmpq_clear(bm);
        fmpq_clear(t);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
