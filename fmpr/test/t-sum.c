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

    printf("sum....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000000; iter++)
    {
        long i, len, prec, bits, expbits, res1, res2;
        fmpr_t s1, s2, s3, err, err_bound;
        fmpr_struct terms[20];
        fmpr_rnd_t rnd;

        len = n_randint(state, 20);
        bits = 2 + n_randint(state, 1000);
        prec = 2 + n_randint(state, 1000);
        expbits = n_randint(state, 14);

        fmpr_init(s1);
        fmpr_init(s2);
        fmpr_init(s3);
        fmpr_init(err);
        fmpr_init(err_bound);
        for (i = 0; i < len; i++)
        {
            fmpr_init(terms + i);
            fmpr_randtest_special(terms + i, state, bits, expbits);
        }

        switch (n_randint(state, 4))
        {
            case 0: rnd = FMPR_RND_DOWN; break;
            case 1: rnd = FMPR_RND_UP; break;
            case 2: rnd = FMPR_RND_FLOOR; break;
            default: rnd = FMPR_RND_CEIL; break;
        }

        res1 = fmpr_sum(s1, terms, len, prec, rnd);

        fmpr_zero(s2);
        for (i = 0; i < len; i++)
            fmpr_add(s2, s2, terms + i, FMPR_PREC_EXACT, FMPR_RND_DOWN);
        res2 = fmpr_set_round(s3, s2, prec, rnd);

        if (!fmpr_equal(s1, s3) || res1 != res2 ||
            !fmpr_check_ulp(s1, res1, prec) || !fmpr_check_ulp(s3, res2, prec))
        {
            printf("FAIL (%ld)\n\n", iter);
            printf("prec = %ld\n\n", prec);
            for (i = 0; i < len; i++)
            {
                printf("terms[%ld] = ", i); fmpr_print(terms + i); printf("\n\n");
            }
            printf("s1 = "); fmpr_print(s1); printf("\n\n");
            printf("s2 = "); fmpr_print(s2); printf("\n\n");
            printf("s3 = "); fmpr_print(s3); printf("\n\n");
            printf("res1 = %ld, res2 = %ld\n\n", res1, res2);
            abort();
        }

        fmpr_sub(err, s1, s2, FMPR_PREC_EXACT, FMPR_RND_DOWN);
        fmpr_abs(err, err);
        fmpr_set_error_result(err_bound, s1, res1);

        if (fmpr_cmp(err, err_bound) > 0)
        {
            printf("FAIL (error bound)!\n");
            printf("prec = %ld\n\n", prec);
            for (i = 0; i < len; i++)
            {
                printf("terms[%ld] = ", i); fmpr_print(terms + i); printf("\n\n");
            }
            printf("s1 = "); fmpr_print(s1); printf("\n\n");
            printf("s2 = "); fmpr_print(s2); printf("\n\n");
            printf("s3 = "); fmpr_print(s3); printf("\n\n");
            printf("error: "); fmpr_print(err); printf("\n\n");
            printf("error bound: "); fmpr_print(err_bound); printf("\n\n");
            printf("res1 = %ld, res2 = %ld\n\n", res1, res2);
            abort();
        }

        fmpr_clear(s1);
        fmpr_clear(s2);
        fmpr_clear(s3);
        fmpr_clear(err);
        fmpr_clear(err_bound);
        for (i = 0; i < len; i++)
            fmpr_clear(terms + i);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

