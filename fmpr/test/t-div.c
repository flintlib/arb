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

    printf("div....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        long bits, r1, r2;
        fmpr_t x, y, z, w;
        mpfr_t X, Y, Z;

        bits = 2 + n_randint(state, 200);

        fmpr_init(x);
        fmpr_init(y);
        fmpr_init(z);
        fmpr_init(w);

        mpfr_init2(X, bits + 100);
        mpfr_init2(Y, bits + 100);
        mpfr_init2(Z, bits);

        fmpr_randtest_special(x, state, bits + n_randint(state, 100), 10);
        fmpr_randtest_special(y, state, bits + n_randint(state, 100), 10);
        fmpr_randtest_special(z, state, bits + n_randint(state, 100), 10);

        fmpr_get_mpfr(X, x, MPFR_RNDN);
        fmpr_get_mpfr(Y, y, MPFR_RNDN);

        switch (n_randint(state, 4))
        {
            case 0:
                r1 = mpfr_div(Z, X, Y, MPFR_RNDZ);
                r2 = fmpr_div(z, x, y, bits, FMPR_RND_DOWN);
                break;
            case 1:
                r1 = mpfr_div(Z, X, Y, MPFR_RNDA);
                r2 = fmpr_div(z, x, y, bits, FMPR_RND_UP);
                break;
            case 2:
                r1 = mpfr_div(Z, X, Y, MPFR_RNDD);
                r2 = fmpr_div(z, x, y, bits, FMPR_RND_FLOOR);
                break;
            default:
                r1 = mpfr_div(Z, X, Y, MPFR_RNDU);
                r2 = fmpr_div(z, x, y, bits, FMPR_RND_CEIL);
                break;
        }

        /* we use slightly different semantics for special values (?) */
        if (mpfr_zero_p(Y))
            mpfr_set_nan(Z);
        if (mpfr_inf_p(X) && mpfr_zero_p(Y))
            mpfr_set_nan(Z);

        fmpr_set_mpfr(w, Z);

        if (!fmpr_equal(z, w) || ((r1 == 0) != (r2 == FMPR_RESULT_EXACT))
            || !fmpr_check_ulp(z, r2, bits))
        {
            printf("FAIL\n\n");
            printf("bits = %ld\n", bits);
            printf("r1 = %ld, r2 = %ld\n", r1, r2);
            printf("x = "); fmpr_print(x); printf("\n\n");
            printf("y = "); fmpr_print(y); printf("\n\n");
            printf("z = "); fmpr_print(z); printf("\n\n");
            printf("w = "); fmpr_print(w); printf("\n\n");
            abort();
        }

        fmpr_clear(x);
        fmpr_clear(y);
        fmpr_clear(z);
        fmpr_clear(w);

        mpfr_clear(X);
        mpfr_clear(Y);
        mpfr_clear(Z);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
