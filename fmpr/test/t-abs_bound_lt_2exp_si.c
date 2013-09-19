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

    printf("abs_bound_lt_2exp_si....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        fmpr_t x;
        fmpz_t b;
        long c;

        fmpr_init(x);
        fmpz_init(b);

        fmpr_randtest_special(x, state, 2 + n_randint(state, 1000), 100);
        fmpr_abs_bound_lt_2exp_fmpz(b, x);
        c = fmpr_abs_bound_lt_2exp_si(x);

        if (c == -FMPR_PREC_EXACT)
        {
            if (!(fmpr_is_zero(x) || fmpz_cmp_si(b, -FMPR_PREC_EXACT) <= 0))
            {
                printf("FAIL (small/zero)\n\n");
                abort();
            }
        }
        else if (c == FMPR_PREC_EXACT)
        {
            if (!(fmpr_is_inf(x) || fmpr_is_nan(x) ||
                fmpz_cmp_si(b, FMPR_PREC_EXACT) >= 0))
            {
                printf("FAIL (large/inf/nan)\n\n");
                abort();
            }
        }
        else
        {
            if (fmpz_cmp_si(b, c) != 0)
            {
                printf("FAIL (normal)\n\n");
                printf("x = "); fmpr_print(x); printf("\n\n");
                printf("b = "); fmpz_print(b); printf("\n\n");
                printf("c = %ld\n\n", c);
                abort();
            }
        }

        fmpr_clear(x);
        fmpz_clear(b);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

