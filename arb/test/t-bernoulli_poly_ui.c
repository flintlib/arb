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

#include <stdio.h>
#include "arb.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("bernoulli_poly_ui....");
    fflush(stdout);

    flint_randinit(state);

    /* test multiplication theorem */
    for (iter = 0; iter < 1000; iter++)
    {
        arb_t x, t, res1, res2;
        ulong n, m, k;
        slong prec;

        n = n_randint(state, 50);
        m = 1 + n_randint(state, 5);
        prec = 2 + n_randint(state, 200);

        arb_init(x);
        arb_init(t);
        arb_init(res1);
        arb_init(res2);

        arb_randtest(x, state, 2 + n_randint(state, 200), 20);
        arb_randtest(res1, state, 2 + n_randint(state, 200), 20);

        arb_mul_ui(t, x, m, prec);
        arb_bernoulli_poly_ui(res1, n, t, prec);

        arb_zero(res2);
        for (k = 0; k < m; k++)
        {
            arb_set_ui(t, k);
            arb_div_ui(t, t, m, prec);
            arb_add(t, t, x, prec);
            arb_bernoulli_poly_ui(t, n, t, prec);
            arb_add(res2, res2, t, prec);
        }

        if (n > 0)
        {
            arb_ui_pow_ui(t, m, n - 1, prec);
            arb_mul(res2, res2, t, prec);
        }
        else
        {
            arb_div_ui(res2, res2, m, prec);
        }

        if (!arb_overlaps(res1, res2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("n = %wu, m = %wu\n\n", n, m);
            flint_printf("x = "); arb_printd(x, 15); flint_printf("\n\n");
            flint_printf("res1 = "); arb_printd(res1, 15); flint_printf("\n\n");
            flint_printf("res2 = "); arb_printd(res2, 15); flint_printf("\n\n");
            abort();
        }

        arb_clear(x);
        arb_clear(t);
        arb_clear(res1);
        arb_clear(res2);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

