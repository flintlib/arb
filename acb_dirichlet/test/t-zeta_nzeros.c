/*
    Copyright (C) 2019 D.H.J Polymath

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

    flint_printf("zeta_nzeros....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 500 * arb_test_multiplier(); iter++)
    {
        arb_t a, b, c;
        slong prec1, prec2;

        prec1 = 2 + n_randint(state, 100);
        prec2 = 2 + n_randint(state, 100);

        arb_init(a);
        arb_init(b);
        arb_init(c);

        arb_randtest_precise(a, state, 1 + n_randint(state, 1000), 3);

        acb_dirichlet_zeta_nzeros(b, a, prec1);

        if (n_randint(state, 2))
        {
            acb_dirichlet_zeta_nzeros(c, a, prec2);
        }
        else  /* test aliasing */
        {
            arb_set(c, a);
            acb_dirichlet_zeta_nzeros(c, c, prec2);
        }

        if (!arb_overlaps(b, c))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
    }

    /* specific examples */
    {
        arb_t t, u, a, b;
        slong i, prec = 53;
        slong exact_int_in[6] = {0, 1, 14, 15, 56, 57};
        slong exact_int_out[6] = {0, 0, 0, 1, 11, 12};

        arb_init(t);
        arb_init(u);
        arb_init(a);
        arb_init(b);

        for (i = 0; i < 6; i++)
        {
            arb_set_si(t, exact_int_in[i]);
            acb_dirichlet_zeta_nzeros(a, t, prec);
            if (!arb_equal_si(a, exact_int_out[i]))
            {
                flint_printf("FAIL: exact small integer\n\n");
                flint_printf("t = "); arb_print(t); flint_printf("\n\n");
                flint_printf("a = "); arb_print(a); flint_printf("\n\n");
                flint_abort();
            }
        }

        arb_set_ui(t, 50);
        mag_one(arb_radref(t));
        arb_set_ui(b, 19);
        mag_one(arb_radref(b));
        arb_mul_2exp_si(b, b, -1);
        acb_dirichlet_zeta_nzeros(a, t, prec);
        if (!arb_equal(a, b))
        {
            flint_printf("FAIL: interval of small integers\n\n");
            flint_printf("t = "); arb_print(t); flint_printf("\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_abort();
        }

        arb_set_str(t, "14.134725141734", prec);
        acb_dirichlet_zeta_nzeros(a, t, prec);
        if (!arb_is_zero(a))
        {
            flint_printf("FAIL: example near first zero\n\n");
            flint_printf("t = "); arb_print(t); flint_printf("\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_abort();
        }

        arb_set_str(t, "14.1347251417346", prec);
        acb_dirichlet_zeta_nzeros(a, t, prec);
        if (!arb_is_zero(a))
        {
            flint_printf("FAIL: example near first zero\n\n");
            flint_printf("t = "); arb_print(t); flint_printf("\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_abort();
        }

        arb_set_str(t, "14.134725141735", prec);
        acb_dirichlet_zeta_nzeros(a, t, prec);
        if (!arb_is_one(a))
        {
            flint_printf("FAIL: example near first zero\n\n");
            flint_printf("t = "); arb_print(t); flint_printf("\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_abort();
        }

        arb_set_ui(t, 10000000);
        acb_dirichlet_zeta_nzeros(a, t, prec);
        if (!arb_is_finite(a) || !arb_contains_si(a, 21136125))
        {
            flint_printf("FAIL: example\n\n");
            flint_printf("t = "); arb_print(t); flint_printf("\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_abort();
        }

        arb_set_str(t, "545439823.215", prec);
        acb_dirichlet_zeta_nzeros(a, t, prec);
        if (!arb_is_finite(a) || !arb_contains_si(a, 1500000001))
        {
            flint_printf("FAIL: example\n\n");
            flint_printf("t = "); arb_print(t); flint_printf("\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(t);
        arb_clear(a);
        arb_clear(b);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
