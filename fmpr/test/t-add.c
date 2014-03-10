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

    Copyright (C) 2012, 2013 Fredrik Johansson

******************************************************************************/

#include "fmpr.h"


int main()
{
    long iter;
    flint_rand_t state;

    printf("add....");
    fflush(stdout);

    flint_randinit(state);

    /* test exact addition: (x + y) + z == x + (y + z) */
    for (iter = 0; iter < 100000; iter++)
    {
        long bits, res1, res2, res3, res4;
        fmpr_t x, y, z, t, u;

        bits = 2 + n_randint(state, 200);

        fmpr_init(x);
        fmpr_init(y);
        fmpr_init(z);
        fmpr_init(t);
        fmpr_init(u);

        fmpr_randtest_special(x, state, bits, 10);
        fmpr_randtest_special(y, state, bits, 10);
        fmpr_randtest_special(z, state, bits, 10);
        fmpr_randtest_special(t, state, bits, 10);
        fmpr_randtest_special(u, state, bits, 10);

        res1 = fmpr_add(t, x, y, FMPR_PREC_EXACT, FMPR_RND_DOWN);
        res2 = fmpr_add(t, t, z, FMPR_PREC_EXACT, FMPR_RND_DOWN);
        res3 = fmpr_add(u, y, z, FMPR_PREC_EXACT, FMPR_RND_DOWN);
        res4 = fmpr_add(u, x, u, FMPR_PREC_EXACT, FMPR_RND_DOWN);

        if (!fmpr_equal(t, u) ||
            res1 != FMPR_RESULT_EXACT || res2 != FMPR_RESULT_EXACT ||
            res3 != FMPR_RESULT_EXACT || res4 != FMPR_RESULT_EXACT)
        {
            printf("FAIL\n\n");
            printf("bits = %ld\n", bits);
            printf("x = "); fmpr_print(x); printf("\n\n");
            printf("y = "); fmpr_print(y); printf("\n\n");
            printf("z = "); fmpr_print(z); printf("\n\n");
            printf("t = "); fmpr_print(t); printf("\n\n");
            printf("u = "); fmpr_print(u); printf("\n\n");
            abort();
        }

        fmpr_clear(x);
        fmpr_clear(y);
        fmpr_clear(z);
        fmpr_clear(t);
        fmpr_clear(u);
    }

    /* compare with add_naive */
    for (iter = 0; iter < 1000000; iter++)
    {
        long prec, ret1, ret2;
        fmpr_t x, y, z, w;
        fmpr_rnd_t rnd;

        fmpr_init(x);
        fmpr_init(y);
        fmpr_init(z);
        fmpr_init(w);

        fmpr_randtest_special(x, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 200));
        fmpr_randtest_special(y, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 200));
        fmpr_randtest_special(z, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 200));

        prec = 2 + n_randint(state, 1000);
        if (n_randint(state, 10) == 0 &&
                fmpz_bits(fmpr_expref(x)) < 10 &&
                fmpz_bits(fmpr_expref(y)) < 10)
            prec = FMPR_PREC_EXACT;

        switch (n_randint(state, 4))
        {
            case 0: rnd = FMPR_RND_DOWN; break;
            case 1: rnd = FMPR_RND_UP; break;
            case 2: rnd = FMPR_RND_FLOOR; break;
            default: rnd = FMPR_RND_CEIL; break;
        }

        ret1 = fmpr_add(z, x, y, prec, rnd);
        ret2 = fmpr_add_naive(w, x, y, prec, rnd);

        if (!fmpr_equal(z, w) || ret1 != ret2 ||
            !fmpr_check_ulp(z, ret1, prec) || !fmpr_check_ulp(w, ret2, prec))
        {
            printf("FAIL\n\n");
            printf("iter %ld\n", iter);
            printf("prec = %ld\n", prec);
            printf("x = "); fmpr_print(x); printf("\n\n");
            printf("y = "); fmpr_print(y); printf("\n\n");
            printf("z = "); fmpr_print(z); printf("\n\n");
            printf("w = "); fmpr_print(w); printf("\n\n");
            printf("ret1 = %ld, ret2 = %ld\n", ret1, ret2);
            abort();
        }

        fmpr_clear(x);
        fmpr_clear(y);
        fmpr_clear(z);
        fmpr_clear(w);
    }

    /* compare rounding with mpfr */
    for (iter = 0; iter < 1000000; iter++)
    {
        long bits, res;
        int mpfr_res;
        fmpr_t x, y, z, w;
        mpfr_t X, Y, Z;

        bits = 2 + n_randint(state, 500);

        fmpr_init(x);
        fmpr_init(y);
        fmpr_init(z);
        fmpr_init(w);

        mpfr_init2(X, bits + 500);
        mpfr_init2(Y, bits + 500);
        mpfr_init2(Z, bits);

        fmpr_randtest_special(x, state, bits + n_randint(state, 500), 10);
        fmpr_randtest_special(y, state, bits + n_randint(state, 500), 10);
        fmpr_randtest_special(z, state, bits, 10);

        fmpr_get_mpfr(X, x, MPFR_RNDN);
        fmpr_get_mpfr(Y, y, MPFR_RNDN);

        switch (n_randint(state, 4))
        {
            case 0:
                mpfr_res = mpfr_add(Z, X, Y, MPFR_RNDZ);
                res = fmpr_add(z, x, y, bits, FMPR_RND_DOWN);
                break;
            case 1:
                mpfr_res = mpfr_add(Z, X, Y, MPFR_RNDA);
                res = fmpr_add(z, x, y, bits, FMPR_RND_UP);
                break;
            case 2:
                mpfr_res = mpfr_add(Z, X, Y, MPFR_RNDD);
                res = fmpr_add(z, x, y, bits, FMPR_RND_FLOOR);
                break;
            default:
                mpfr_res = mpfr_add(Z, X, Y, MPFR_RNDU);
                res = fmpr_add(z, x, y, bits, FMPR_RND_CEIL);
                break;
        }

        fmpr_set_mpfr(w, Z);

        if (!fmpr_equal(z, w) || (res == FMPR_RESULT_EXACT) != (mpfr_res == 0)
            || !fmpr_check_ulp(z, res, bits))
        {
            printf("FAIL\n\n");
            printf("iter %ld\n", iter);
            printf("bits = %ld\n", bits);
            printf("x = "); fmpr_print(x); printf("\n\n");
            printf("y = "); fmpr_print(y); printf("\n\n");
            printf("z = "); fmpr_print(z); printf("\n\n");
            printf("w = "); fmpr_print(w); printf("\n\n");
            printf("returned: %ld, %d\n", res, mpfr_res);
            abort();
        }

        /* check error bound */
        if (!fmpr_is_nan(z) && !fmpr_is_inf(z))
        {
            fmpr_t z_exact, error_bound, true_error;

            fmpr_init(z_exact);
            fmpr_init(error_bound);
            fmpr_init(true_error);

            fmpr_set_error_result(error_bound, z, res);
            fmpr_add(z_exact, x, y, FMPR_PREC_EXACT, FMPR_RND_DOWN);

            fmpr_sub(true_error, z, z_exact, FMPR_PREC_EXACT, FMPR_RND_DOWN);
            fmpr_abs(true_error, true_error);

            if (fmpr_is_zero(error_bound) != fmpr_is_zero(true_error) ||
                fmpr_cmp(true_error, error_bound) > 0)
            {
                printf("FAIL: error bound\n\n");
                printf("bits = %ld\n", bits);
                printf("x = "); fmpr_print(x); printf("\n\n");
                printf("y = "); fmpr_print(y); printf("\n\n");
                printf("z = "); fmpr_print(z); printf("\n\n");
                printf("z_exact = "); fmpr_print(z_exact); printf("\n\n");
                printf("true_error = "); fmpr_print(true_error); printf("\n\n");
                printf("error_bound = "); fmpr_print(error_bound); printf("\n\n");
                abort();
            }

            fmpr_clear(z_exact);
            fmpr_clear(error_bound);
            fmpr_clear(true_error);
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
