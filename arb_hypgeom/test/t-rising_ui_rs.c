/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_hypgeom.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("rising_ui_rs....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * arb_test_multiplier(); iter++)
    {
        arb_t x, xk, y, ya, yb, yayb;
        ulong k, n, m1, m2, m3;
        slong prec;

        prec = 2 + n_randint(state, 200);
        k = n_randint(state, 10);
        n = n_randint(state, 50);
        m1 = n_randint(state, 2) ? 0 : 1 + n_randint(state, FLINT_MAX(n + k, 1));
        m2 = n_randint(state, 2) ? 0 : 1 + n_randint(state, FLINT_MAX(k, 1));
        m3 = n_randint(state, 2) ? 0 : 1 + n_randint(state, FLINT_MAX(n, 1));

        arb_init(x);
        arb_init(xk);
        arb_init(y);
        arb_init(ya);
        arb_init(yb);
        arb_init(yayb);

        arb_randtest(x, state, prec, 10);
        arb_add_ui(xk, x, k, prec);

        arb_hypgeom_rising_ui_rs(y, x, n + k, m1, prec);
        arb_hypgeom_rising_ui_rs(ya, x, k, m2, prec);
        arb_hypgeom_rising_ui_rs(yb, xk, n, m3, prec);
        arb_mul(yayb, ya, yb, prec);

        if (!arb_overlaps(y, yayb))
        {
            flint_printf("FAIL\n\n");
            flint_printf("k = %wu, n = %wu, m1 = %wu, m2 = %wu, m3 = %wu\n\n", k, n, m1, m2, m3);
            flint_printf("x = "); arb_printn(x, 100, 0); flint_printf("\n\n");
            flint_printf("y = "); arb_printn(y, 100, 0); flint_printf("\n\n");
            flint_printf("ya = "); arb_printn(ya, 100, 0); flint_printf("\n\n");
            flint_printf("yb = "); arb_printn(yb, 100, 0); flint_printf("\n\n");
            flint_printf("yayb = "); arb_printn(yayb, 100, 0); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x);
        arb_clear(xk);
        arb_clear(y);
        arb_clear(ya);
        arb_clear(yb);
        arb_clear(yayb);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
