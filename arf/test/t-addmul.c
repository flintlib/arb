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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "arf.h"

int
arf_addmul_naive(arf_t z, const arf_t x, const arf_t y, long prec, arf_rnd_t rnd)
{
    arf_t t;
    int inexact;

    arf_init(t);
    arf_mul(t, x, y, ARF_PREC_EXACT, ARF_RND_DOWN);

    inexact = arf_add(z, z, t, prec, rnd);

    arf_clear(t);

    return inexact;
}

int main()
{
    long iter, iter2;
    flint_rand_t state;

    printf("addmul....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        arf_t x, y, z, v;
        long prec, r1, r2;
        arf_rnd_t rnd;

        arf_init(x);
        arf_init(y);
        arf_init(z);
        arf_init(v);

        for (iter2 = 0; iter2 < 100; iter2++)
        {
            arf_randtest_special(x, state, 2000, 100);
            arf_randtest_special(y, state, 2000, 100);
            arf_randtest_special(z, state, 2000, 100);
            arf_set(v, z);

            prec = 2 + n_randint(state, 2000);

            if (n_randint(state, 10) == 0 &&
                fmpz_bits(ARF_EXPREF(x)) < 10 &&
                fmpz_bits(ARF_EXPREF(y)) < 10 &&
                fmpz_bits(ARF_EXPREF(z)) < 10)
            {
                prec = ARF_PREC_EXACT;
            }

            switch (n_randint(state, 4))
            {
                case 0:  rnd = ARF_RND_DOWN; break;
                case 1:  rnd = ARF_RND_UP; break;
                case 2:  rnd = ARF_RND_FLOOR; break;
                default: rnd = ARF_RND_CEIL; break;
            }

            switch (n_randint(state, 5))
            {
            case 0:
                r1 = arf_addmul(z, x, y, prec, rnd);
                r2 = arf_addmul_naive(v, x, y, prec, rnd);
                if (!arf_equal(z, v) || r1 != r2)
                {
                    printf("FAIL!\n");
                    printf("prec = %ld, rnd = %d\n\n", prec, rnd);
                    printf("x = "); arf_print(x); printf("\n\n");
                    printf("y = "); arf_print(y); printf("\n\n");
                    printf("z = "); arf_print(z); printf("\n\n");
                    printf("v = "); arf_print(v); printf("\n\n");
                    printf("r1 = %ld, r2 = %ld\n", r1, r2);
                    abort();
                }
                break;

            case 1:
                r1 = arf_addmul(z, x, x, prec, rnd);
                r2 = arf_addmul_naive(v, x, x, prec, rnd);
                if (!arf_equal(z, v) || r1 != r2)
                {
                    printf("FAIL (aliasing 1)!\n");
                    printf("prec = %ld, rnd = %d\n\n", prec, rnd);
                    printf("x = "); arf_print(x); printf("\n\n");
                    printf("z = "); arf_print(z); printf("\n\n");
                    printf("v = "); arf_print(v); printf("\n\n");
                    printf("r1 = %ld, r2 = %ld\n", r1, r2);
                    abort();
                }
                break;

            case 2:
                r2 = arf_addmul_naive(v, v, v, prec, rnd);
                r1 = arf_addmul(z, z, z, prec, rnd);
                if (!arf_equal(v, z) || r1 != r2)
                {
                    printf("FAIL (aliasing 2)!\n");
                    printf("prec = %ld, rnd = %d\n\n", prec, rnd);
                    printf("v = "); arf_print(v); printf("\n\n");
                    printf("z = "); arf_print(z); printf("\n\n");
                    printf("r1 = %ld, r2 = %ld\n", r1, r2);
                    abort();
                }
                break;

            case 3:
                r2 = arf_addmul_naive(v, v, y, prec, rnd);
                r1 = arf_addmul(z, z, y, prec, rnd);
                if (!arf_equal(v, z) || r1 != r2)
                {
                    printf("FAIL (aliasing 3)!\n");
                    printf("prec = %ld, rnd = %d\n\n", prec, rnd);
                    printf("y = "); arf_print(y); printf("\n\n");
                    printf("v = "); arf_print(v); printf("\n\n");
                    printf("z = "); arf_print(z); printf("\n\n");
                    printf("r1 = %ld, r2 = %ld\n", r1, r2);
                    abort();
                }
                break;

            default:
                r2 = arf_addmul_naive(v, x, v, prec, rnd);
                r1 = arf_addmul(z, x, z, prec, rnd);
                if (!arf_equal(z, v) || r1 != r2)
                {
                    printf("FAIL (aliasing 4)!\n");
                    printf("prec = %ld, rnd = %d\n\n", prec, rnd);
                    printf("x = "); arf_print(x); printf("\n\n");
                    printf("v = "); arf_print(v); printf("\n\n");
                    printf("z = "); arf_print(z); printf("\n\n");
                    printf("r1 = %ld, r2 = %ld\n", r1, r2);
                    abort();
                }
                break;
            }
        }

        arf_clear(x);
        arf_clear(y);
        arf_clear(z);
        arf_clear(v);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
