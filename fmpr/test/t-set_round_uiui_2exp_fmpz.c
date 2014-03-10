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

    printf("set_round_uiui_2exp_fmpz....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        long prec, ret1, ret2;
        fmpz_t man, exp;
        fmpr_t x, y;
        mp_limb_t hi, lo;
        fmpr_rnd_t rnd;
        int negative;

        fmpz_init(man);
        fmpz_init(exp);
        fmpr_init(x);
        fmpr_init(y);

        prec = 2 + n_randint(state, 1000);
        if (n_randint(state, 10) == 0)
            prec = FMPR_PREC_EXACT;

        negative = n_randint(state, 2);
        hi = n_randtest(state);
        lo = n_randtest(state);
        fmpz_randtest(exp, state, 1 + n_randint(state, 100));
        if (negative)
            fmpz_neg_uiui(man, hi, lo);
        else
            fmpz_set_uiui(man, hi, lo);

        switch (n_randint(state, 4))
        {
            case 0: rnd = FMPR_RND_DOWN; break;
            case 1: rnd = FMPR_RND_UP; break;
            case 2: rnd = FMPR_RND_FLOOR; break;
            default: rnd = FMPR_RND_CEIL; break;
        }

        ret1 = fmpr_set_round_uiui_2exp_fmpz(x, hi, lo, exp, negative, prec, rnd);

        fmpr_set_fmpz_2exp(y, man, exp);
        ret2 = fmpr_set_round(y, y, prec, rnd);

        if (!fmpr_equal(x, y) || ret1 != ret2 ||
            !fmpr_check_ulp(x, ret1, prec) || !fmpr_check_ulp(y, ret2, prec))
        {
            printf("FAIL\n\n");
            printf("prec: %ld\n", prec);
            printf("hi = %lu, lo = %lu\n", hi, lo);
            printf("man = "); fmpz_print(man); printf("\n\n");
            printf("exp = "); fmpz_print(exp); printf("\n\n");
            printf("x = "); fmpr_print(x); printf("\n\n");
            printf("y = "); fmpr_print(y); printf("\n\n");
            printf("ret1 = %ld, ret2 = %ld\n\n", ret1, ret2);
            abort();
        }

        fmpz_clear(man);
        fmpz_clear(exp);
        fmpr_clear(x);
        fmpr_clear(y);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

