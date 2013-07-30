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

#include "fmpz_holonomic.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("forward_fmpz_mat....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        fmpz_holonomic_t op;
        fmpz_mat_t M1, M2, M3, L, R;
        fmpz_t Q1, Q2, Q3;
        long start, n1, n2, r;

        fmpz_holonomic_init(op);
        fmpz_holonomic_randtest(op, state, 5, 5, 10);

        r = fmpz_holonomic_order(op);
        start = n_randint(state, 20);
        n1 = n_randint(state, 10);
        n2 = n_randint(state, 10);

        fmpz_mat_init(M1, r, r);
        fmpz_mat_init(M2, r, r);
        fmpz_mat_init(M3, r, r);
        fmpz_mat_init(L, r, r);
        fmpz_mat_init(R, r, r);
        fmpz_init(Q1);
        fmpz_init(Q2);
        fmpz_init(Q3);

        fmpz_holonomic_forward_fmpz_mat(M1, Q1, op, start, n1);
        fmpz_holonomic_forward_fmpz_mat(M2, Q2, op, start + n1, n2);
        fmpz_holonomic_forward_fmpz_mat(M3, Q3, op, start, n1 + n2);

        /* test Q3 * M2 * M1 = M3 * Q1 * Q2 */
        fmpz_mat_mul(L, M2, M1);
        fmpz_mat_scalar_mul_fmpz(L, L, Q3);

        fmpz_mat_scalar_mul_fmpz(R, M3, Q1);
        fmpz_mat_scalar_mul_fmpz(R, R, Q2);

        if (!fmpz_mat_equal(L, R))
        {
            printf("FAIL\n");

            printf("start = %ld, n1 = %ld, n2 = %ld\n", start, n1, n2);
            fmpz_holonomic_print(op, "n", "Sn"); printf("\n\n");

            fmpz_mat_print_pretty(M1); printf("\n\n");
            fmpz_mat_print_pretty(M2); printf("\n\n");
            fmpz_mat_print_pretty(M3); printf("\n\n");

            abort();
        }

        fmpz_mat_clear(M1);
        fmpz_mat_clear(M2);
        fmpz_mat_clear(M3);
        fmpz_mat_clear(L);
        fmpz_mat_clear(R);
        fmpz_clear(Q1);
        fmpz_clear(Q2);
        fmpz_clear(Q3);
        fmpz_holonomic_clear(op);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

