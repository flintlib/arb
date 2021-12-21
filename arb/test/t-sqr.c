/*
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

int
mag_close(const mag_t am, const mag_t bm)
{
    arf_t t, a, b;
    int res1, res2;

    arf_init(t);
    arf_init(a);
    arf_init(b);

    arf_set_mag(a, am);
    arf_set_mag(b, bm);

    arf_mul_ui(t, b, 257, MAG_BITS, ARF_RND_UP);
    arf_mul_2exp_si(t, t, -8);
    res1 = arf_cmp(a, t) <= 0;

    arf_mul_ui(t, a, 257, MAG_BITS, ARF_RND_UP);
    arf_mul_2exp_si(t, t, -8);
    res2 = arf_cmp(b, t) <= 0;

    arf_clear(t);
    arf_clear(a);
    arf_clear(b);

    return res1 && res2;
}

int main()
{
    slong iter, iter2;
    flint_rand_t state;

    flint_printf("sqr....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        arb_t x, z, v;
        slong prec;

        arb_init(x);
        arb_init(z);
        arb_init(v);

        for (iter2 = 0; iter2 < 100; iter2++)
        {
            arb_randtest_special(x, state, n_randint(state,2) ? 2000 : 200, 200);

            prec = 2 + n_randint(state, 2000);

            switch (n_randint(state, 2))
            {
            case 0:
                arb_mul(z, x, x, prec);
                arb_sqr(v, x, prec);

                if (!arf_equal(arb_midref(z), arb_midref(v))
                    || !mag_close(arb_radref(z), arb_radref(v)))
                {
                    flint_printf("FAIL!\n");
                    flint_printf("x = "); arb_print(x); flint_printf("\n\n");
                    flint_printf("z = "); arb_print(z); flint_printf("\n\n");
                    flint_printf("v = "); arb_print(v); flint_printf("\n\n");
                    flint_abort();
                }
                break;

            default:
                arb_mul(z, x, x, prec);
                arb_sqr(x, x, prec);

                if (!arf_equal(arb_midref(x), arb_midref(z))
                    || !mag_close(arb_radref(x), arb_radref(z)))
                {
                    flint_printf("FAIL (aliasing 4)!\n");
                    flint_printf("x = "); arb_print(x); flint_printf("\n\n");
                    flint_printf("v = "); arb_print(v); flint_printf("\n\n");
                    flint_printf("z = "); arb_print(z); flint_printf("\n\n");
                    flint_abort();
                }
                break;
            }
        }

        arb_clear(x);
        arb_clear(z);
        arb_clear(v);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
