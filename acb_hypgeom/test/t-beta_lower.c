/*
    Copyright (C) 2016 Fredrik Johansson

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

    flint_printf("beta_lower....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * arb_test_multiplier(); iter++)
    {
        acb_t a, a1, b, b1, z, w, wa1, wb1, t, u;
        slong prec;
        int regularized;

        acb_init(a);
        acb_init(a1);
        acb_init(b);
        acb_init(b1);
        acb_init(z);
        acb_init(w);
        acb_init(wa1);
        acb_init(wb1);
        acb_init(t);
        acb_init(u);

        regularized = n_randint(state, 2);

        prec = 2 + n_randint(state, 100);

        acb_randtest_param(a, state, 1 + n_randint(state, 100), 1 + n_randint(state, 10));
        acb_randtest_param(b, state, 1 + n_randint(state, 100), 1 + n_randint(state, 10));
        acb_randtest_param(z, state, 1 + n_randint(state, 100), 1 + n_randint(state, 10));

        acb_add_ui(a1, a, 1, prec);
        acb_add_ui(b1, b, 1, prec);

        acb_hypgeom_beta_lower(w, a, b, z, regularized, prec);
        acb_hypgeom_beta_lower(wa1, a1, b, z, regularized, prec);
        acb_hypgeom_beta_lower(wb1, a, b1, z, regularized, prec);

        if (regularized)
        {
            acb_mul(t, wa1, a, prec);
            acb_addmul(t, wb1, b, prec);
            acb_add(u, a, b, prec);
            acb_div(t, t, u, prec);
        }
        else
        {
            acb_add(t, wa1, wb1, prec);
        }

        if (!acb_overlaps(w, t))
        {
            flint_printf("FAIL: contiguous relation\n\n");
            flint_printf("regularized = %d\n\n", regularized);
            flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
            flint_printf("b = "); acb_printd(b, 30); flint_printf("\n\n");
            flint_printf("z = ");  acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("w = "); acb_printd(w, 30); flint_printf("\n\n");
            flint_printf("wa1 = "); acb_printd(wa1, 30); flint_printf("\n\n");
            flint_printf("wb1 = "); acb_printd(wb1, 30); flint_printf("\n\n");
            flint_printf("t = "); acb_printd(t, 30); flint_printf("\n\n");
            flint_abort();
        }

        /* test I(a,b;z) = 1-I(b,a,1-z) */

        if (!regularized)
        {
            acb_add(t, a, b, prec);
            acb_gamma(t, t, prec);
            acb_mul(w, w, t, prec);
            acb_rgamma(t, a, prec);
            acb_mul(w, w, t, prec);
            acb_rgamma(t, b, prec);
            acb_mul(w, w, t, prec);
        }

        acb_sub_ui(t, z, 1, prec);
        acb_neg(t, t);
        acb_hypgeom_beta_lower(t, b, a, t, 1, prec);
        acb_sub_ui(t, t, 1, prec);
        acb_neg(t, t);

        if (!acb_overlaps(w, t))
        {
            flint_printf("FAIL: symmetry\n\n");
            flint_printf("regularized = %d\n\n", regularized);
            flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
            flint_printf("b = "); acb_printd(b, 30); flint_printf("\n\n");
            flint_printf("z = ");  acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("w = "); acb_printd(w, 30); flint_printf("\n\n");
            flint_printf("t = "); acb_printd(t, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(a);
        acb_clear(a1);
        acb_clear(b);
        acb_clear(b1);
        acb_clear(z);
        acb_clear(w);
        acb_clear(wa1);
        acb_clear(wb1);
        acb_clear(t);
        acb_clear(u);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

