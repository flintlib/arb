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
    long iter;
    flint_rand_t state;

    printf("set_round_mpn....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000000; iter++)
    {
        long prec, bits, shift, ret1, ret2;
        fmpz_t man, exp;
        fmpr_t y;
        mpz_t x;
        fmpz_t fx;
        fmpr_rnd_t rnd;

        fmpz_init(man);
        fmpz_init(exp);
        mpz_init(x);
        fmpz_init(fx);
        fmpr_init(y);

        bits = 1 + n_randint(state, 1000);
        prec = 2 + n_randint(state, 1000);
        if (n_randint(state, 10) == 0)
            prec = FMPR_PREC_EXACT;

        fmpz_randtest_not_zero(fx, state, bits);
        fmpz_get_mpz(x, fx);

        switch (n_randint(state, 4))
        {
            case 0: rnd = FMPR_RND_DOWN; break;
            case 1: rnd = FMPR_RND_UP; break;
            case 2: rnd = FMPR_RND_FLOOR; break;
            default: rnd = FMPR_RND_CEIL; break;
        }

        ret1 = _fmpr_set_round_mpn(&shift, man, x->_mp_d,
            FLINT_ABS(x->_mp_size), (x->_mp_size < 0) ? 1 : 0, prec, rnd);
        fmpz_set_si(exp, shift);

        ret2 = fmpr_set_round_fmpz(y, fx, prec, rnd);

        if (!fmpz_equal(fmpr_manref(y), man) ||
            !fmpz_equal(fmpr_expref(y), exp) || ret1 != ret2 ||
            !fmpr_check_ulp(y, ret2, prec))
        {
            printf("FAIL\n\n");
            printf("bits: %ld\n", bits);
            printf("prec: %ld\n", prec);
            printf("x = "); fmpz_print(fx); printf("\n\n");
            printf("man = "); fmpz_print(man); printf("\n\n");
            printf("exp = "); fmpz_print(exp); printf("\n\n");
            printf("y = "); fmpr_print(y); printf("\n\n");
            printf("ret1 = %ld, ret2 = %ld\n\n", ret1, ret2);
            abort();
        }

        fmpz_clear(man);
        fmpz_clear(exp);
        mpz_clear(x);
        fmpz_clear(fx);
        fmpr_clear(y);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
