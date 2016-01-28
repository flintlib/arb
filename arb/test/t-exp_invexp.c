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

    Copyright (C) 2012, 2013 Fredrik Johansson

******************************************************************************/

#include "arb.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("exp_invexp....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        arb_t a, b, c, d;
        slong prec;

        if (iter % 10 == 0)
            prec = 10000;
        else
            prec = 1000;

        prec = 2 + n_randint(state, prec);

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arb_init(d);

        arb_randtest_special(a, state, 1 + n_randint(state, prec), 1 + n_randint(state, 100));
        arb_randtest_special(b, state, 1 + n_randint(state, prec), 1 + n_randint(state, 100));
        arb_randtest_special(c, state, 1 + n_randint(state, prec), 1 + n_randint(state, 100));

        arb_exp_invexp(b, c, a, prec);

        arb_exp(d, a, prec);

        if (!arb_overlaps(b, d))
        {
            flint_printf("FAIL: overlap 1\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_printf("d = "); arb_print(c); flint_printf("\n\n");
            abort();
        }

        arb_neg(d, a);
        arb_exp(d, d, prec);

        if (!arb_overlaps(c, d))
        {
            flint_printf("FAIL: overlap 2\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_printf("d = "); arb_print(c); flint_printf("\n\n");
            abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        arb_clear(d);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

