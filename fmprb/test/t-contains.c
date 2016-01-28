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

    flint_printf("contains....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        fmprb_t a, b;
        fmpq_t am, ar, bm, br, t, u;
        int c1, c2;

        fmprb_init(a);
        fmprb_init(b);

        fmpq_init(am);
        fmpq_init(ar);
        fmpq_init(bm);
        fmpq_init(br);
        fmpq_init(t);
        fmpq_init(u);

        fmprb_randtest(a, state, 1 + n_randint(state, 500), 14);
        fmprb_randtest(b, state, 1 + n_randint(state, 500), 14);

        fmpr_get_fmpq(am, fmprb_midref(a));
        fmpr_get_fmpq(ar, fmprb_radref(a));
        fmpr_get_fmpq(bm, fmprb_midref(b));
        fmpr_get_fmpq(br, fmprb_radref(b));

        c1 = fmprb_contains(a, b);

        fmpq_sub(t, am, ar);
        fmpq_sub(u, bm, br);
        c2 = fmpq_cmp(t, u) <= 0;

        fmpq_add(t, am, ar);
        fmpq_add(u, bm, br);
        c2 = c2 && (fmpq_cmp(t, u) >= 0);

        if (c1 != c2)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "); fmprb_print(a); flint_printf("\n\n");
            flint_printf("b = "); fmprb_print(b); flint_printf("\n\n");
            flint_printf("am = "); fmpq_print(am); flint_printf("\n\n");
            flint_printf("ar = "); fmpq_print(ar); flint_printf("\n\n");
            flint_printf("bm = "); fmpq_print(bm); flint_printf("\n\n");
            flint_printf("br = "); fmpq_print(br); flint_printf("\n\n");
            flint_printf("t = "); fmpq_print(t); flint_printf("\n\n");
            flint_printf("u = "); fmpq_print(u); flint_printf("\n\n");
            flint_printf("c1 = %d, c2 = %d\n\n", c1, c2);
            abort();
        }

        fmprb_clear(a);
        fmprb_clear(b);

        fmpq_clear(am);
        fmpq_clear(ar);
        fmpq_clear(bm);
        fmpq_clear(br);
        fmpq_clear(t);
        fmpq_clear(u);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
