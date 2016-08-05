/*
    Copyright (C) 2014 Fredrik Johansson
    Copyright (C) 2016 Ricky E. Farr

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/


#include <arf.h>

int
arf_fam_naive(arf_t w, const arf_t z, const arf_t x, const arf_t y,
              slong prec, arf_rnd_t rnd)
{
    int inexact;
    arf_t t;
    arf_init(t);
    arf_mul(t, x, y, ARF_PREC_EXACT, ARF_RND_DOWN);
    inexact = arf_add(w, z, t, prec, rnd);
    return inexact;
}

int main(void)
{
    slong iter, iter2;
    slong prec, r1, r2;
    flint_rand_t state;
    arf_t v, x, y, z;
    arf_rnd_t rnd;

    flint_printf("fam....");
    fflush(stdout);

    flint_randinit(state);

    arf_init(v);
    arf_init(x);
    arf_init(y);
    arf_init(z);

    for (iter = 0; iter < 1000 * arb_test_multiplier(); iter++)
    {
        for (iter2 = 0; iter2 < 100; iter2++)
        {
            arf_randtest_special(x, state, 2000, 100);
            arf_randtest_special(y, state, 2000, 100);
            arf_randtest_special(z, state, 2000, 100);

            prec = 2 + n_randint(state, 2000);

            if (n_randint(state, 10) == 0 &&
                fmpz_bits(ARF_EXPREF(x)) < 10 &&
                fmpz_bits(ARF_EXPREF(y)) < 10 &&
                fmpz_bits(ARF_EXPREF(z)) < 10)
            {
                prec = ARF_PREC_EXACT;
            }

            switch (n_randint(state, 5))
            {
                case 0:  rnd = ARF_RND_DOWN; break;
                case 1:  rnd = ARF_RND_UP; break;
                case 2:  rnd = ARF_RND_FLOOR; break;
                case 3:  rnd = ARF_RND_CEIL; break;
                default: rnd = ARF_RND_NEAR; break;
            }

            switch (n_randint(state, 6))
            {
            case 0:
                arf_set(v, z);
                r2 = arf_fam_naive(v, z, x, y, prec, rnd);
                r1 = arf_fam(z, z, x, y, prec, rnd);

                if (!arf_equal(z, v) || r1 != r2)
                {
                    flint_printf("FAIL!\n");
                    flint_printf("prec = %wd, rnd = %d\n\n", prec, rnd);
                    flint_printf("x = "); arf_print(x); flint_printf("\n\n");
                    flint_printf("y = "); arf_print(y); flint_printf("\n\n");
                    flint_printf("z = "); arf_print(z); flint_printf("\n\n");
                    flint_printf("v = "); arf_print(v); flint_printf("\n\n");
                    flint_printf("r1 = %wd, r2 = %wd\n", r1, r2);
                    abort();
                }

                break;

            case 1:
                arf_set(v, z);
                r2 = arf_fam_naive(v, z, x, x, prec, rnd);
                r1 = arf_fam(z, z, x, x, prec, rnd);

                if (!arf_equal(z, v) || r1 != r2)
                {
                    flint_printf("FAIL (aliasing 1)!\n");
                    flint_printf("prec = %wd, rnd = %d\n\n", prec, rnd);
                    flint_printf("x = "); arf_print(x); flint_printf("\n\n");
                    flint_printf("z = "); arf_print(z); flint_printf("\n\n");
                    flint_printf("v = "); arf_print(v); flint_printf("\n\n");
                    flint_printf("r1 = %wd, r2 = %wd\n", r1, r2);
                    abort();
                }

                break;

            case 2:
                arf_set(v, z);
                r2 = arf_fam_naive(v, z, z, z, prec, rnd);
                r1 = arf_fam(z, z, z, z, prec, rnd);
                if (!arf_equal(z, v) || r1 != r2)
                {
                    flint_printf("FAIL (aliasing 2)!\n");
                    flint_printf("prec = %wd, rnd = %d\n\n", prec, rnd);
                    flint_printf("z = "); arf_print(z); flint_printf("\n\n");
                    flint_printf("v = "); arf_print(v); flint_printf("\n\n");
                    flint_printf("r1 = %wd, r2 = %wd\n", r1, r2);
                    abort();
                }

                break;

            case 3:
                arf_set(v, z);
                r2 = arf_fam_naive(v, z, z, y, prec, rnd);
                r1 = arf_fam(z, z, z, y, prec, rnd);

                if (!arf_equal(v, z) || r1 != r2)
                {
                    flint_printf("FAIL (aliasing 3)!\n");
                    flint_printf("prec = %wd, rnd = %d\n\n", prec, rnd);
                    flint_printf("y = "); arf_print(y); flint_printf("\n\n");
                    flint_printf("z = "); arf_print(z); flint_printf("\n\n");
                    flint_printf("v = "); arf_print(v); flint_printf("\n\n");
                    flint_printf("r1 = %wd, r2 = %wd\n", r1, r2);
                    abort();
                }

                break;

            case 4:

                r2 = arf_fam_naive(v, z, z, y, prec, rnd);
                r1 = arf_fam(z, z, z, y, prec, rnd);

                if (!arf_equal(v, z) || r1 != r2)
                {
                    flint_printf("FAIL (v != z)!\n");
                    flint_printf("prec = %wd, rnd = %d\n\n", prec, rnd);
                    flint_printf("y = "); arf_print(y); flint_printf("\n\n");
                    flint_printf("z = "); arf_print(z); flint_printf("\n\n");
                    flint_printf("v = "); arf_print(v); flint_printf("\n\n");
                    flint_printf("r1 = %wd, r2 = %wd\n", r1, r2);
                    abort();
                }

                break;

            default:
                arf_set(v, z);
                r2 = arf_fam_naive(v, z, x, z, prec, rnd);
                r1 = arf_fam(z, z, x, z, prec, rnd);

                if (!arf_equal(z, v) || r1 != r2)
                {
                    flint_printf("FAIL (aliasing 4)!\n");
                    flint_printf("prec = %wd, rnd = %d\n\n", prec, rnd);
                    flint_printf("x = "); arf_print(x); flint_printf("\n\n");
                    flint_printf("v = "); arf_print(v); flint_printf("\n\n");
                    flint_printf("z = "); arf_print(z); flint_printf("\n\n");
                    flint_printf("r1 = %wd, r2 = %wd\n", r1, r2);
                    abort();
                }

                break;
            }
        }
    }

    flint_printf("PASS\n");
    flint_randclear(state);
    arf_clear(v);
    arf_clear(x);
    arf_clear(y);
    arf_clear(z);
    flint_cleanup();

    return EXIT_SUCCESS;
}
