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

    printf("sqrt....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        long bits, res;
        fmpr_t x, z, w;
        mpfr_t X, Z;

        bits = 2 + n_randint(state, 200);

        fmpr_init(x);
        fmpr_init(z);
        fmpr_init(w);

        mpfr_init2(X, 2 * (bits + 100));
        mpfr_init2(Z, bits);

        fmpr_randtest_special(x, state, bits + n_randint(state, 100), 10);
        fmpr_randtest_special(z, state, bits + n_randint(state, 100), 10);

        /* occasionally produce perfect squares */
        if (n_randint(state, 4) == 0)
            fmpr_mul(x, x, x, 2 * bits, FMPR_RND_DOWN);

        fmpr_get_mpfr(X, x, MPFR_RNDN);

        switch (n_randint(state, 4))
        {
            case 0:
                mpfr_sqrt(Z, X, MPFR_RNDZ);
                res = fmpr_sqrt(z, x, bits, FMPR_RND_DOWN);
                break;
            case 1:
                mpfr_sqrt(Z, X, MPFR_RNDA);
                res = fmpr_sqrt(z, x, bits, FMPR_RND_UP);
                break;
            case 2:
                mpfr_sqrt(Z, X, MPFR_RNDD);
                res = fmpr_sqrt(z, x, bits, FMPR_RND_FLOOR);
                break;
            default:
                mpfr_sqrt(Z, X, MPFR_RNDU);
                res = fmpr_sqrt(z, x, bits, FMPR_RND_CEIL);
                break;
        }

        fmpr_set_mpfr(w, Z);

        if (!fmpr_equal(z, w) || !fmpr_check_ulp(z, res, bits))
        {
            printf("FAIL\n\n");
            printf("bits = %ld\n", bits);
            printf("x = "); fmpr_print(x); printf("\n\n");
            printf("z = "); fmpr_print(z); printf("\n\n");
            printf("w = "); fmpr_print(w); printf("\n\n");
            abort();
        }

        fmpr_clear(x);
        fmpr_clear(z);
        fmpr_clear(w);

        mpfr_clear(X);
        mpfr_clear(Z);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
