/*
    Copyright (C) 2012, 2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arf.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("get_mpfr....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000 * arb_test_multiplier(); iter++)
    {
        slong bits;
        arf_t x, z;
        mpfr_t y;

        bits = 2 + n_randint(state, 200);

        arf_init(x);
        arf_init(z);
        mpfr_init2(y, bits);

        arf_randtest_special(x, state, bits, 10);
        arf_get_mpfr(y, x, MPFR_RNDN);
        arf_set_mpfr(z, y);

        if (!arf_equal(x, z))
        {
            flint_printf("FAIL\n\n");
            flint_printf("x = "); arf_print(x); flint_printf("\n\n");
            flint_printf("z = "); arf_print(z); flint_printf("\n\n");
            abort();
        }

        arf_clear(x);
        arf_clear(z);
        mpfr_clear(y);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

