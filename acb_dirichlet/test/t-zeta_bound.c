/*
    Copyright (C) 2016 Fredrik Johansson

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

    flint_printf("zeta_bound....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * arb_test_multiplier(); iter++)
    {
        acb_t s, sb, z;
        arb_t x, y;
        slong prec, i;
        mag_t zbound;

        acb_init(s);
        acb_init(sb);
        acb_init(z);
        arb_init(x);
        arb_init(y);
        mag_init(zbound);

        acb_randtest(s, state, 2 + n_randint(state, 100), 4);
        if (n_randint(state, 2))
            arb_set_d(acb_realref(s), -0.5 + n_randint(state, 32) / 16.0);
        mag_zero(arb_radref(acb_realref(s)));
        mag_zero(arb_radref(acb_imagref(s)));

        /* create a larger interval containing s */
        acb_set(sb, s);
        for (i = 0; i < 2; i++)
        {
            arb_randtest(x, state, 2 + n_randint(state, 100), 2);
            arb_add(x, x, acb_realref(s), 2 * MAG_BITS);
            arb_intersection(acb_realref(sb), acb_realref(sb), x, 2 * MAG_BITS);

            arb_randtest(x, state, 2 + n_randint(state, 100), 2);
            arb_add(x, x, acb_imagref(s), 2 * MAG_BITS);
            arb_intersection(acb_imagref(sb), acb_imagref(sb), x, 2 * MAG_BITS);
        }

        acb_dirichlet_zeta_bound(zbound, sb);

        if (!mag_is_inf(zbound))
        {
            for (prec = 64; ; prec *= 2)
            {
                acb_zeta(z, s, prec);
                if (acb_rel_accuracy_bits(z) > MAG_BITS)
                    break;
            }

            acb_abs(x, z, 2 * MAG_BITS);

            arb_zero(y);
            arf_set_mag(arb_midref(y), zbound);

            if (!arb_le(x, y))
            {
                flint_printf("FAIL: bound\n\n");
                flint_printf("iter = %wd\n", iter);
                flint_printf("s = "); acb_printn(s, 50, 0); flint_printf("\n\n");
                flint_printf("sb = "); acb_printn(sb, 50, 0); flint_printf("\n\n");
                flint_printf("z = "); acb_printn(z, 50, 0); flint_printf("\n\n");
                flint_printf("zbound = "); mag_printd(zbound, 10); flint_printf("\n\n");
                flint_printf("x = "); arb_printn(x, 50, 0); flint_printf("\n\n");
                flint_printf("y = "); arb_printn(y, 50, 0); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_clear(s);
        acb_clear(sb);
        acb_clear(z);
        arb_clear(x);
        arb_clear(y);
        mag_clear(zbound);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

