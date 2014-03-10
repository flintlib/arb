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

#include "fmpr.h"

int main()
{
    long iter, iter2;
    flint_rand_t state;

    printf("mul....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        fmpr_t x, y, z, v;
        long prec, r1, r2;
        fmpr_rnd_t rnd;

        fmpr_init(x);
        fmpr_init(y);
        fmpr_init(z);
        fmpr_init(v);

        for (iter2 = 0; iter2 < 100; iter2++)
        {
            fmpr_randtest_special(x, state, 2000, 200);
            fmpr_randtest_special(y, state, 2000, 200);
            prec = 2 + n_randint(state, 2000);

            switch (n_randint(state, 4))
            {
                case 0:  rnd = FMPR_RND_DOWN; break;
                case 1:  rnd = FMPR_RND_UP; break;
                case 2:  rnd = FMPR_RND_FLOOR; break;
                default: rnd = FMPR_RND_CEIL; break;
            }

            switch (n_randint(state, 5))
            {
            case 0:
                r1 = fmpr_mul(z, x, y, prec, rnd);
                r2 = fmpr_mul_naive(v, x, y, prec, rnd);
                if (!fmpr_equal(z, v) || r1 != r2 || !fmpr_check_ulp(z, r1, prec))
                {
                    printf("FAIL!\n");
                    printf("x = "); fmpr_print(x); printf("\n\n");
                    printf("y = "); fmpr_print(y); printf("\n\n");
                    printf("z = "); fmpr_print(z); printf("\n\n");
                    printf("v = "); fmpr_print(v); printf("\n\n");
                    printf("r1 = %ld, r2 = %ld\n", r1, r2);
                    abort();
                }
                break;

            case 1:
                r1 = fmpr_mul(z, x, x, prec, rnd);
                r2 = fmpr_mul_naive(v, x, x, prec, rnd);
                if (!fmpr_equal(z, v) || r1 != r2 || !fmpr_check_ulp(z, r1, prec))
                {
                    printf("FAIL (aliasing 1)!\n");
                    printf("x = "); fmpr_print(x); printf("\n\n");
                    printf("z = "); fmpr_print(z); printf("\n\n");
                    printf("v = "); fmpr_print(v); printf("\n\n");
                    printf("r1 = %ld, r2 = %ld\n", r1, r2);
                    abort();
                }
                break;

            case 2:
                r2 = fmpr_mul_naive(v, x, x, prec, rnd);
                r1 = fmpr_mul(x, x, x, prec, rnd);
                if (!fmpr_equal(v, x) || r1 != r2 || !fmpr_check_ulp(x, r1, prec))
                {
                    printf("FAIL (aliasing 2)!\n");
                    printf("x = "); fmpr_print(x); printf("\n\n");
                    printf("z = "); fmpr_print(z); printf("\n\n");
                    printf("v = "); fmpr_print(v); printf("\n\n");
                    printf("r1 = %ld, r2 = %ld\n", r1, r2);
                    abort();
                }
                break;

            case 3:
                r2 = fmpr_mul_naive(v, x, y, prec, rnd);
                r1 = fmpr_mul(x, x, y, prec, rnd);
                if (!fmpr_equal(x, v) || r1 != r2 || !fmpr_check_ulp(x, r1, prec))
                {
                    printf("FAIL (aliasing 3)!\n");
                    printf("x = "); fmpr_print(x); printf("\n\n");
                    printf("y = "); fmpr_print(y); printf("\n\n");
                    printf("v = "); fmpr_print(v); printf("\n\n");
                    printf("r1 = %ld, r2 = %ld\n", r1, r2);
                    abort();
                }
                break;

            default:
                r2 = fmpr_mul_naive(v, x, y, prec, rnd);
                r1 = fmpr_mul(x, y, x, prec, rnd);
                if (!fmpr_equal(x, v) || r1 != r2 || !fmpr_check_ulp(x, r1, prec))
                {
                    printf("FAIL (aliasing 4)!\n");
                    printf("x = "); fmpr_print(x); printf("\n\n");
                    printf("y = "); fmpr_print(y); printf("\n\n");
                    printf("v = "); fmpr_print(v); printf("\n\n");
                    printf("r1 = %ld, r2 = %ld\n", r1, r2);
                    abort();
                }
                break;
            }
        }

        fmpr_clear(x);
        fmpr_clear(y);
        fmpr_clear(z);
        fmpr_clear(v);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
