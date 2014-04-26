/*=============================================================================

    This file is part of ARB.

    ARB is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    ARB is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ARB; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2014 Fredrik Johansson

******************************************************************************/

#include "arf.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("set_round....");
    fflush(stdout);

    flint_randinit(state);

    {
        arf_t x, y, z;

        arf_init(x);
        arf_init(y);
        arf_init(z);

        for (iter = 0; iter < 100000; iter++)
        {
            long bits1, bits2;
            int ret1, ret2;
            mpfr_t g1, g2;
            fmpz_t e;
            arf_rnd_t rnd;

            bits1 = 1 + n_randint(state, 1000);
            bits2 = 2 + n_randint(state, 1000);

            if (n_randint(state, 100) == 0)
                bits2 = ARF_PREC_EXACT;

            switch (n_randint(state, 4))
            {
                case 0: rnd = ARF_RND_DOWN; break;
                case 1: rnd = ARF_RND_UP; break;
                case 2: rnd = ARF_RND_FLOOR; break;
                default: rnd = ARF_RND_CEIL; break;
            }

            fmpz_init(e);
            mpfr_init2(g1, FLINT_MAX(2, bits1));
            mpfr_init2(g2, FLINT_MIN(bits2, 10000));

            if (n_randint(state, 100) == 0)
            {
                arf_clear(x); arf_clear(y); arf_clear(z);
                arf_init(x); arf_init(y); arf_init(z);
            }

            /* dirty output variables */
            if (n_randint(state, 2))
            {
                arf_randtest_special(y, state, 1 + n_randint(state, 1000),
                    1 + n_randint(state, 100));
                arf_randtest_special(z, state, 1 + n_randint(state, 1000),
                    1 + n_randint(state, 100));
            }

            arf_randtest_special(x, state, bits1, 1 + n_randint(state, 10));
            arf_get_mpfr(g1, x, MPFR_RNDD); /* exact */

            /* test large exponents */
            if (n_randint(state, 4) == 0)
                fmpz_randtest(e, state, 1 + n_randint(state, 100));

            if (!arf_is_special(x))
                fmpz_add(ARF_EXPREF(x), ARF_EXPREF(x), e);

            ret1 = arf_set_round(y, x, bits2, rnd);
            ret2 = mpfr_set(g2, g1, arf_rnd_to_mpfr(rnd));
            arf_set_mpfr(z, g2);

            if (!arf_is_special(y))
                fmpz_sub(ARF_EXPREF(y), ARF_EXPREF(y), e);

            if (!arf_equal(y, z) || ((ret1 == ARF_RESULT_EXACT) != (ret2 == 0)))
            {
                printf("FAIL\n\n");
                printf("bits1: %ld\n", bits1);
                printf("bits2: %ld\n", bits2);
                printf("x = "); arf_print(x); printf("\n\n");
                printf("y = "); arf_print(y); printf("\n\n");
                printf("z = "); arf_print(z); printf("\n\n");
                printf("ret1 = %d, ret2 = %d\n\n", ret1, ret2);
                abort();
            }

            if (!arf_is_special(x))
                fmpz_add(ARF_EXPREF(x), ARF_EXPREF(x), e);

            ret1 = arf_set_round(y, x, bits2, rnd);
            arf_set(z, x);
            ret2 = arf_set_round(z, z, bits2, rnd);

            if (!arf_equal(y, z) || ret1 != ret2)
            {
                printf("FAIL (aliasing)\n\n");
                printf("bits1: %ld\n", bits1);
                printf("bits2: %ld\n", bits2);
                printf("x = "); arf_print(x); printf("\n\n");
                printf("y = "); arf_print(y); printf("\n\n");
                printf("z = "); arf_print(z); printf("\n\n");
                printf("ret1 = %d, ret2 = %d\n\n", ret1, ret2);
                abort();
            }

            mpfr_clear(g1);
            mpfr_clear(g2);
            fmpz_clear(e);
        }

        arf_clear(x);
        arf_clear(y);
        arf_clear(z);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

