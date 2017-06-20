/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("sinc_pi....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * arb_test_multiplier(); iter++)
    {
        arb_t x, y, z;
        slong prec;

        arb_init(x);
        arb_init(y);
        arb_init(z);

        prec = 2 + n_randint(state, 200);

        arb_randtest(x, state, 1 + n_randint(state, 200), 1 + n_randint(state, 100));
        arb_randtest(y, state, 1 + n_randint(state, 200), 1 + n_randint(state, 100));

        if (n_randint(state, 2))
        {
            arb_set(y, x);
            arb_sinc_pi(y, y, prec);
        }
        else
            arb_sinc_pi(y, x, prec);

        arb_const_pi(z, prec);
        arb_mul(z, z, x, prec);
        arb_sinc(z, z, prec);

        if (!arb_overlaps(y, z))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("x = "); arb_printd(x, 30); flint_printf("\n\n");
            flint_printf("y = "); arb_printd(y, 30); flint_printf("\n\n");
            flint_printf("z = "); arb_printd(z, 30); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x);
        arb_clear(y);
        arb_clear(z);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
