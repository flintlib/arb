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

#include "fmpr.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("ulp....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        fmpr_t x, ulp, a, b;
        long prec;

        fmpr_init(x);
        fmpr_init(ulp);
        fmpr_init(a);
        fmpr_init(b);

        prec = 2 + n_randint(state, 1000);
        fmpr_randtest_not_zero(x, state, 1 + n_randint(state, 1000),
            1 + n_randint(state, 100));

        fmpr_ulp(ulp, x, prec);

        fmpr_abs(a, x);
        fmpr_mul_2exp_si(a, a, -prec);
        fmpr_abs(b, x);
        fmpr_mul_2exp_si(b, b, -prec+1);

        if (!((fmpr_cmp(a, ulp) < 0) && (fmpr_cmp(ulp, b) <= 0)))
        {
            printf("FAIL!\n");
            printf("prec = %ld\n", prec);
            printf("x = "); fmpr_print(x); printf("\n");
            printf("ulp = "); fmpr_print(ulp); printf("\n");
            printf("a = "); fmpr_print(a); printf("\n");
            printf("b = "); fmpr_print(b); printf("\n");
            abort();
        }

        fmpr_clear(x);
        fmpr_clear(ulp);
        fmpr_clear(a);
        fmpr_clear(b);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

