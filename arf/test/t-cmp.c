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

    Copyright (C) 2012, 2014 Fredrik Johansson

******************************************************************************/

#include "arf.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("cmp....");
    fflush(stdout);

    flint_randinit(state);

    /* compare with fmpz */
    {
        arf_t x, y;
        fmpz_t X, Y;

        arf_init(x);
        arf_init(y);
        fmpz_init(X);
        fmpz_init(Y);

        for (iter = 0; iter < 100000; iter++)
        {
            int cmp1, cmp2;

            fmpz_randtest(X, state, 1 + n_randint(state, 1000));

            switch (n_randint(state, 8))
            {
                case 0:
                    fmpz_neg(Y, X);
                    break;
                case 1:
                    fmpz_set(Y, X);
                    break;
                default:
                    fmpz_randtest(Y, state, 1 + n_randint(state, 1000));
            }

            arf_set_fmpz(x, X);
            arf_set_fmpz(y, Y);

            cmp1 = arf_cmp(x, y);

            cmp2 = fmpz_cmp(X, Y);
            cmp2 = (cmp2 > 0) - (cmp2 < 0);

            if (cmp1 != cmp2)
            {
                printf("FAIL\n\n");
                printf("x = "); arf_debug(x); printf("\n\n");
                printf("y = "); arf_debug(y); printf("\n\n");
                printf("X = "); fmpz_print(X); printf("\n\n");
                printf("Y = "); fmpz_print(Y); printf("\n\n");
                printf("cmp1 = %d, cmp2 = %d\n\n", cmp1, cmp2);
                abort();
            }
        }

        arf_clear(x);
        arf_clear(y);
        fmpz_clear(X);
        fmpz_clear(Y);
    }

    /* compare with mpfr */
    for (iter = 0; iter < 100000; iter++)
    {
        long bits;
        arf_t x, y;
        mpfr_t X, Y;
        int cmp1, cmp2;

        bits = 2 + n_randint(state, 200);

        arf_init(x);
        arf_init(y);

        mpfr_init2(X, bits);
        mpfr_init2(Y, bits);

        arf_randtest_special(x, state, bits, 10);
        arf_randtest_special(y, state, bits, 10);

        arf_get_mpfr(X, x, MPFR_RNDN);
        arf_get_mpfr(Y, y, MPFR_RNDN);

        cmp1 = arf_cmp(x, y);
        cmp2 = mpfr_cmp(X, Y);

        if (cmp1 != cmp2)
        {
            printf("FAIL\n\n");
            printf("x = "); arf_print(x); printf("\n\n");
            printf("y = "); arf_print(y); printf("\n\n");
            printf("cmp1 = %d, cmp2 = %d\n\n", cmp1, cmp2);
            abort();
        }

        arf_clear(x);
        arf_clear(y);

        mpfr_clear(X);
        mpfr_clear(Y);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

