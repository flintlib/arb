/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_hypgeom.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("rising_ui_rs....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * arb_test_multiplier(); iter++)
    {
        acb_t x, xk, y, ya, yb, yayb;
        ulong k, n, m1, m2, m3;
        slong prec;

        prec = 2 + n_randint(state, 200);
        k = n_randint(state, 10);
        n = n_randint(state, 50);
        m1 = n_randint(state, 2) ? 0 : 1 + n_randint(state, FLINT_MAX(n + k, 1));
        m2 = n_randint(state, 2) ? 0 : 1 + n_randint(state, FLINT_MAX(k, 1));
        m3 = n_randint(state, 2) ? 0 : 1 + n_randint(state, FLINT_MAX(n, 1));

        acb_init(x);
        acb_init(xk);
        acb_init(y);
        acb_init(ya);
        acb_init(yb);
        acb_init(yayb);

        acb_randtest(x, state, prec, 10);
        acb_add_ui(xk, x, k, prec);

        acb_hypgeom_rising_ui_rs(y, x, n + k, m1, prec);
        acb_hypgeom_rising_ui_rs(ya, x, k, m2, prec);
        acb_hypgeom_rising_ui_rs(yb, xk, n, m3, prec);
        acb_mul(yayb, ya, yb, prec);

        if (!acb_overlaps(y, yayb))
        {
            flint_printf("FAIL\n\n");
            flint_printf("k = %wu, n = %wu, m1 = %wu, m2 = %wu, m3 = %wu\n\n", k, n, m1, m2, m3);
            flint_printf("x = "); acb_printn(x, 100, 0); flint_printf("\n\n");
            flint_printf("y = "); acb_printn(y, 100, 0); flint_printf("\n\n");
            flint_printf("ya = "); acb_printn(ya, 100, 0); flint_printf("\n\n");
            flint_printf("yb = "); acb_printn(yb, 100, 0); flint_printf("\n\n");
            flint_printf("yayb = "); acb_printn(yayb, 100, 0); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(x);
        acb_clear(xk);
        acb_clear(y);
        acb_clear(ya);
        acb_clear(yb);
        acb_clear(yayb);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
