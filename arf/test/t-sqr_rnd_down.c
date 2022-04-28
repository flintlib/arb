/*
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arf.h"

int main()
{
    slong iter, iter2;
    flint_rand_t state;

    flint_printf("sqr_rnd_down....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        arf_t x, z, v;
        slong prec, r1, r2;
        arf_rnd_t rnd;

        arf_init(x);
        arf_init(z);
        arf_init(v);

        for (iter2 = 0; iter2 < 100; iter2++)
        {
            arf_randtest_special(x, state, 2000, 100);
            prec = 2 + n_randint(state, 2000);

            if (n_randint(state, 50) == 0)
                prec = ARF_PREC_EXACT;

            r1 = arf_mul_rnd_down(z, x, x, prec);
            r2 = arf_sqr_rnd_down(v, x, prec);
            if (!arf_equal(z, v) || r1 != r2)
            {
                flint_printf("FAIL!\n");
                flint_printf("prec = %wd, rnd = %d\n\n", prec, rnd);
                flint_printf("x = "); arf_print(x); flint_printf("\n\n");
                flint_printf("z = "); arf_print(z); flint_printf("\n\n");
                flint_printf("v = "); arf_print(v); flint_printf("\n\n");
                flint_printf("r1 = %wd, r2 = %wd\n", r1, r2);
                flint_abort();
            }
        }

        arf_clear(x);
        arf_clear(z);
        arf_clear(v);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
