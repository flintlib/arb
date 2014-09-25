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

    Copyright (C) 2014 Fredrik Johansson

******************************************************************************/

#include "arb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("get_interval_mpfr....");
    fflush(stdout);
    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        arb_t x, y;
        mpfr_t aa, bb;

        arb_init(x);
        arb_init(y);
        mpfr_init2(aa, 2 + n_randint(state, 200));
        mpfr_init2(bb, 2 + n_randint(state, 200));

        arb_randtest_special(x, state, 200, 10);
        arb_get_interval_mpfr(aa, bb, x);
        arb_set_interval_mpfr(y, aa, bb, 2 + n_randint(state, 200));

        if (!arb_contains(y, x))
        {
            printf("FAIL:\n\n");
            printf("x = "); arb_print(x); printf("\n\n");
            printf("aa = "); mpfr_printf("%.50Rg", aa); printf("\n\n");
            printf("bb = "); mpfr_printf("%.50Rg", bb); printf("\n\n");
            printf("y = "); arb_print(y); printf("\n\n");
            abort();
        }

        arb_clear(x);
        arb_clear(y);
        mpfr_clear(aa);
        mpfr_clear(bb);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

