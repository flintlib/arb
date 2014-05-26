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

#include "arf.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("sum....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000000; iter++)
    {
        long i, len, prec, bits, expbits;
        int res1, res2;
        arf_t s1, s2, s3, err;
        mag_t err_bound;
        arf_struct terms[20];
        arf_rnd_t rnd;

        len = n_randint(state, 20);
        bits = 2 + n_randint(state, 1000);
        prec = 2 + n_randint(state, 1000);
        expbits = n_randint(state, 14);

        arf_init(s1);
        arf_init(s2);
        arf_init(s3);
        arf_init(err);
        mag_init(err_bound);

        for (i = 0; i < len; i++)
        {
            arf_init(terms + i);
            arf_randtest_special(terms + i, state, bits, expbits);
        }

        switch (n_randint(state, 4))
        {
            case 0: rnd = ARF_RND_DOWN; break;
            case 1: rnd = ARF_RND_UP; break;
            case 2: rnd = ARF_RND_FLOOR; break;
            default: rnd = ARF_RND_CEIL; break;
        }

        res1 = arf_sum(s1, terms, len, prec, rnd);

        arf_zero(s2);
        for (i = 0; i < len; i++)
            arf_add(s2, s2, terms + i, ARF_PREC_EXACT, ARF_RND_DOWN);
        res2 = arf_set_round(s3, s2, prec, rnd);

        if (!arf_equal(s1, s3) || res1 != res2)
        {
            printf("FAIL (%ld)\n\n", iter);
            printf("prec = %ld\n\n", prec);
            for (i = 0; i < len; i++)
            {
                printf("terms[%ld] = ", i); arf_print(terms + i); printf("\n\n");
            }
            printf("s1 = "); arf_print(s1); printf("\n\n");
            printf("s2 = "); arf_print(s2); printf("\n\n");
            printf("s3 = "); arf_print(s3); printf("\n\n");
            printf("res1 = %d, res2 = %d\n\n", res1, res2);
            abort();
        }

        arf_sub(err, s1, s2, ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_abs(err, err);

        if (res1)
            arf_mag_set_ulp(err_bound, s1, prec);
        else
            mag_zero(err_bound);

        if (arf_cmpabs_mag(err, err_bound) > 0)
        {
            printf("FAIL (error bound)!\n");
            printf("prec = %ld\n\n", prec);
            for (i = 0; i < len; i++)
            {
                printf("terms[%ld] = ", i); arf_print(terms + i); printf("\n\n");
            }
            printf("s1 = "); arf_print(s1); printf("\n\n");
            printf("s2 = "); arf_print(s2); printf("\n\n");
            printf("s3 = "); arf_print(s3); printf("\n\n");
            printf("error: "); arf_print(err); printf("\n\n");
            printf("error bound: "); mag_print(err_bound); printf("\n\n");
            printf("res1 = %d, res2 = %d\n\n", res1, res2);
            abort();
        }

        arf_clear(s1);
        arf_clear(s2);
        arf_clear(s3);
        arf_clear(err);
        mag_clear(err_bound);

        for (i = 0; i < len; i++)
            arf_clear(terms + i);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

