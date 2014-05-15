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

#include "acb_mat.h"


int main()
{
    long iter;
    flint_rand_t state;

    printf("mul....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        long m, n, k, qbits1, qbits2, rbits1, rbits2, rbits3;
        fmpq_mat_t A, B, C;
        acb_mat_t a, b, c, d;

        qbits1 = 2 + n_randint(state, 200);
        qbits2 = 2 + n_randint(state, 200);
        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);
        rbits3 = 2 + n_randint(state, 200);

        m = n_randint(state, 10);
        n = n_randint(state, 10);
        k = n_randint(state, 10);

        fmpq_mat_init(A, m, n);
        fmpq_mat_init(B, n, k);
        fmpq_mat_init(C, m, k);

        acb_mat_init(a, m, n);
        acb_mat_init(b, n, k);
        acb_mat_init(c, m, k);
        acb_mat_init(d, m, k);

        fmpq_mat_randtest(A, state, qbits1);
        fmpq_mat_randtest(B, state, qbits2);
        fmpq_mat_mul(C, A, B);

        acb_mat_set_fmpq_mat(a, A, rbits1);
        acb_mat_set_fmpq_mat(b, B, rbits2);
        acb_mat_mul(c, a, b, rbits3);

        if (!acb_mat_contains_fmpq_mat(c, C))
        {
            printf("FAIL\n\n");
            printf("m = %ld, n = %ld, k = %ld, bits3 = %ld\n", m, n, k, rbits3);

            printf("A = "); fmpq_mat_print(A); printf("\n\n");
            printf("B = "); fmpq_mat_print(B); printf("\n\n");
            printf("C = "); fmpq_mat_print(C); printf("\n\n");

            printf("a = "); acb_mat_printd(a, 15); printf("\n\n");
            printf("b = "); acb_mat_printd(b, 15); printf("\n\n");
            printf("c = "); acb_mat_printd(c, 15); printf("\n\n");

            abort();
        }

        /* test aliasing with a */
        if (acb_mat_nrows(a) == acb_mat_nrows(c) &&
            acb_mat_ncols(a) == acb_mat_ncols(c))
        {
            acb_mat_set(d, a);
            acb_mat_mul(d, d, b, rbits3);
            if (!acb_mat_equal(d, c))
            {
                printf("FAIL (aliasing 1)\n\n");
                abort();
            }
        }

        /* test aliasing with b */
        if (acb_mat_nrows(b) == acb_mat_nrows(c) &&
            acb_mat_ncols(b) == acb_mat_ncols(c))
        {
            acb_mat_set(d, b);
            acb_mat_mul(d, a, d, rbits3);
            if (!acb_mat_equal(d, c))
            {
                printf("FAIL (aliasing 2)\n\n");
                abort();
            }
        }

        fmpq_mat_clear(A);
        fmpq_mat_clear(B);
        fmpq_mat_clear(C);

        acb_mat_clear(a);
        acb_mat_clear(b);
        acb_mat_clear(c);
        acb_mat_clear(d);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
