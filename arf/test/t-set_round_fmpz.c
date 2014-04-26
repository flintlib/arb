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

    printf("set_round_fmpz....");
    fflush(stdout);

    flint_randinit(state);

    {
        arf_t x, y;

        arf_init(x);
        arf_init(y);

        for (iter = 0; iter < 100000; iter++)
        {
            long bits1, bits2;
            int ret1, ret2;
            fmpz_t a;
            mpz_t b;
            mpfr_t g;
            arf_rnd_t rnd;

            bits1 = 1 + n_randint(state, 1000);
            bits2 = 2 + n_randint(state, 1000);

            if (n_randint(state, 100) == 0)
                bits2 = ARF_PREC_EXACT;

            fmpz_init(a);
            mpz_init(b);
            mpfr_init2(g, FLINT_MIN(bits2, 10000));

            if (n_randint(state, 100) == 0)
            {
                arf_clear(x);
                arf_clear(y);
                arf_init(x);
                arf_init(y);
            }

            /* dirty output variables */
            if (n_randint(state, 2))
            {
                arf_randtest_special(x, state, 1 + n_randint(state, 1000),
                    1 + n_randint(state, 100));
                arf_randtest_special(y, state, 1 + n_randint(state, 1000),
                    1 + n_randint(state, 100));
            }

            fmpz_randtest(a, state, bits1);
            fmpz_get_mpz(b, a);

            switch (n_randint(state, 4))
            {
                case 0: rnd = ARF_RND_DOWN; break;
                case 1: rnd = ARF_RND_UP; break;
                case 2: rnd = ARF_RND_FLOOR; break;
                default: rnd = ARF_RND_CEIL; break;
            }

            ret1 = arf_set_round_fmpz(x, a, bits2, rnd);
            ret2 = mpfr_set_z(g, b, arf_rnd_to_mpfr(rnd));
            arf_set_mpfr(y, g);
            arf_equal(x, y);

            if (!arf_equal(x, y) || ((ret1 == ARF_RESULT_EXACT) != (ret2 == 0)))
            {
                printf("FAIL\n\n");
                printf("bits1: %ld\n", bits1);
                printf("bits2: %ld\n", bits2);
                printf("a = "); fmpz_print(a); printf("\n\n");
                printf("x = "); arf_print(x); printf("\n\n");
                printf("y = "); arf_print(y); printf("\n\n");
                printf("ret1 = %d, ret2 = %d\n\n", ret1, ret2);
                abort();
            }

            fmpz_clear(a);
            mpz_clear(b);
            mpfr_clear(g);
        }

        arf_clear(x);
        arf_clear(y);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

