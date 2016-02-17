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

    Copyright (C) 2016 Fredrik Johansson

******************************************************************************/

#include "arb.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("sgn....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        arb_t a, b;
        int result;

        arb_init(a);
        arb_init(b);

        arb_randtest_special(a, state, 1 + n_randint(state, 200), 10);
        arb_randtest_special(b, state, 1 + n_randint(state, 200), 10);
        arb_sgn(b, a);

        result = 1;
        if (arb_contains_zero(a))
            result = result & arb_contains_si(b, 0);
        if (arb_contains_positive(a))
            result = result & arb_contains_si(b, 1);
        if (arb_contains_negative(a))
            result = result & arb_contains_si(b, -1);

        if (!result)
        {
            flint_printf("FAIL\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            abort();
        }

        arb_clear(a);
        arb_clear(b);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
