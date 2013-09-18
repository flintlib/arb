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

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include "fmpr.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("abs_bound_lt_2exp_fmpz....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        fmpr_t x, y;
        fmpz_t b;
        int cmp1, cmp2;

        fmpr_init(x);
        fmpr_init(y);
        fmpz_init(b);

        fmpr_randtest_not_zero(x, state, 2 + n_randint(state, 1000), 100);
        fmpr_abs_bound_lt_2exp_fmpz(b, x);

        fmpr_one(y);
        fmpr_mul_2exp_fmpz(y, y, b);

        cmp1 = (fmpr_cmpabs(x, y) < 0);

        fmpr_mul_2exp_si(y, y, -1);

        cmp2 = (fmpr_cmpabs(y, x) <= 0);

        fmpr_mul_2exp_si(y, y, 1);

        if (!cmp1 || !cmp2)
        {
            printf("FAIL\n\n");
            printf("x = "); fmpr_print(x); printf("\n\n");
            printf("y = "); fmpr_print(y); printf("\n\n");
            printf("b = "); fmpz_print(b); printf("\n\n");
            printf("cmp1 = %d, cmp2 = %d\n\n", cmp1, cmp2);
            abort();
        }

        fmpr_clear(x);
        fmpr_clear(y);
        fmpz_clear(b);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

