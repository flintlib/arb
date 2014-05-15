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
    long iter, iter2;
    flint_rand_t state;

    printf("complex_sqr....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        arf_t e1, f1, e2, f2, a, b;
        long prec, r1, r2;
        arf_rnd_t rnd;

        arf_init(a);
        arf_init(b);
        arf_init(e1);
        arf_init(f1);
        arf_init(e2);
        arf_init(f2);

        for (iter2 = 0; iter2 < 100; iter2++)
        {
            arf_randtest_special(a, state, 3000, 100);
            arf_randtest_special(b, state, 3000, 100);
            prec = 2 + n_randint(state, 3000);

            switch (n_randint(state, 4))
            {
                case 0:  rnd = ARF_RND_DOWN; break;
                case 1:  rnd = ARF_RND_UP; break;
                case 2:  rnd = ARF_RND_FLOOR; break;
                default: rnd = ARF_RND_CEIL; break;
            }

            switch (n_randint(state, 2))
            {
            case 0:
                r1 = arf_complex_sqr(e1, f1, a, b, prec, rnd);
                r2 = arf_complex_mul_fallback(e2, f2, a, b, a, b, prec, rnd);
                if (!arf_equal(e1, e2) || !arf_equal(f1, f2) || r1 != r2)
                {
                    printf("FAIL!\n");
                    printf("prec = %ld, rnd = %d\n\n", prec, rnd);
                    printf("a = "); arf_print(a); printf("\n\n");
                    printf("b = "); arf_print(b); printf("\n\n");
                    printf("e1 = "); arf_print(e1); printf("\n\n");
                    printf("f1 = "); arf_print(f1); printf("\n\n");
                    printf("e2 = "); arf_print(e2); printf("\n\n");
                    printf("f2 = "); arf_print(f2); printf("\n\n");
                    printf("r1 = %ld, r2 = %ld\n", r1, r2);
                    abort();
                }
                break;

            default:
                r1 = arf_complex_mul_fallback(e1, f1, a, b, a, b, prec, rnd);
                r2 = arf_complex_sqr(a, b, a, b, prec, rnd);
                if (!arf_equal(e1, a) || !arf_equal(f1, b) || r1 != r2)
                {
                    printf("FAIL! (aliasing)\n");
                    printf("prec = %ld, rnd = %d\n\n", prec, rnd);
                    printf("a = "); arf_print(a); printf("\n\n");
                    printf("b = "); arf_print(b); printf("\n\n");
                    printf("e1 = "); arf_print(e1); printf("\n\n");
                    printf("f1 = "); arf_print(f1); printf("\n\n");
                    printf("r1 = %ld, r2 = %ld\n", r1, r2);
                    abort();
                }
                break;
            }
        }

        arf_clear(a);
        arf_clear(b);
        arf_clear(e1);
        arf_clear(f1);
        arf_clear(e2);
        arf_clear(f2);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
