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

    printf("set_interval_mpfr....");
    fflush(stdout);
    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        arb_t x;
        arf_t a, b;
        mpfr_t aa, bb;

        arb_init(x);
        arf_init(a);
        arf_init(b);
        mpfr_init2(aa, 200);
        mpfr_init2(bb, 200);

        arf_randtest_special(a, state, 200, 10);
        arf_randtest_special(b, state, 200, 10);
        if (arf_cmp(a, b) > 0)
            arf_swap(a, b);

        arf_get_mpfr(aa, a, MPFR_RNDD);
        arf_get_mpfr(bb, b, MPFR_RNDU);

        arb_set_interval_mpfr(x, aa, bb, 2 + n_randint(state, 200));

        if (!arb_contains_arf(x, a) || !arb_contains_arf(x, b))
        {
            printf("FAIL:\n\n");
            printf("x = "); arb_print(x); printf("\n\n");
            printf("a = "); arf_print(a); printf("\n\n");
            printf("b = "); arf_print(b); printf("\n\n");
            abort();
        }

        arb_clear(x);
        arf_clear(a);
        arf_clear(b);
        mpfr_clear(aa);
        mpfr_clear(bb);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

