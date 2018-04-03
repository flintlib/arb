/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("stieltjes....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 250 * arb_test_multiplier(); iter++)
    {
        acb_t b1, b2, a;
        fmpz_t n;
        slong prec1, prec2, acc1, acc2;

        fmpz_init(n);
        acb_init(b1);
        acb_init(b2);
        acb_init(a);

        if (n_randint(state, 10) == 0)
        {
            acb_one(a);
            fmpz_randtest(n, state, 1 + n_randint(state, 100));
            fmpz_abs(n, n);
            prec1 = 2 + n_randint(state, 100);
            prec2 = 2 + n_randint(state, 100);
        }
        else if (n_randint(state, 7) != 0)
        {
            acb_one(a);
            fmpz_randtest(n, state, 2 + n_randint(state, 12));
            fmpz_abs(n, n);
            prec1 = 2 + n_randint(state, 400);
            prec2 = 2 + n_randint(state, 400);
        }
        else
        {
            acb_randtest_precise(a, state, 400, 2);
            fmpz_randtest(n, state, 7);
            fmpz_abs(n, n);
            prec1 = 2 + n_randint(state, 400);
            prec2 = 2 + n_randint(state, 400);
        }

        acb_dirichlet_stieltjes(b1, n, a, prec1);
        acb_dirichlet_stieltjes(b2, n, a, prec2);

        if (!acb_overlaps(b1, b2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("n = "); fmpz_print(n); flint_printf("\n\n");
            flint_printf("b1 = "); acb_printn(b1, 1000, 0); flint_printf("\n\n");
            flint_printf("b2 = "); acb_printn(b2, 1000, 0); flint_printf("\n\n");
            flint_abort();
        }

        acc1 = acb_rel_accuracy_bits(b1);
        acc2 = acb_rel_accuracy_bits(b2);

        if (acb_is_one(a) && (acc1 < prec1 - 10 || acc2 < prec2 - 10))
        {
            flint_printf("FAIL: poor accuracy\n\n");
            flint_printf("prec1 = %wd, acc1 = %wd\n", prec1, acc1);
            flint_printf("prec2 = %wd, acc2 = %wd\n", prec2, acc2);
            flint_printf("n = "); fmpz_print(n); flint_printf("\n\n");
            flint_printf("b1 = "); acb_printn(b1, 500, 0); flint_printf("\n\n");
            flint_printf("b2 = "); acb_printn(b2, 500, 0); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(b1);
        acb_clear(b2);
        acb_clear(a);
        fmpz_clear(n);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

