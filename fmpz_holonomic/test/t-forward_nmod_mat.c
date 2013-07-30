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

    printf("forward_nmod_mat....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        fmpz_holonomic_t op;
        fmpz_mat_t M1;
        fmpz_t Q1;
        nmod_mat_t M2, L, R;
        mp_limb_t Q2;
        long start, n, r;
        mp_limb_t p;

        fmpz_holonomic_init(op);

        fmpz_holonomic_randtest(op, state, 4, 4, 10);

        r = fmpz_holonomic_order(op);
        start = n_randint(state, 10);
        n = n_randint(state, 100);
        p = n_randtest_prime(state, 0);

        fmpz_mat_init(M1, r, r);
        fmpz_init(Q1);
        nmod_mat_init(M2, r, r, p);
        nmod_mat_init(L, r, r, p);
        nmod_mat_init(R, r, r, p);

        fmpz_holonomic_forward_fmpz_mat(M1, Q1, op, start, n);
        fmpz_holonomic_forward_nmod_mat(M2, &Q2, op, start, n);

        fmpz_mat_get_nmod_mat(L, M1);
        nmod_mat_scalar_mul(L, L, Q2);

        nmod_mat_scalar_mul(R, M2, fmpz_fdiv_ui(Q1, p));

        /* check Q2 * M1 = Q1 * M2 */

        if (!nmod_mat_equal(L, R))
        {
            printf("FAIL\n");

            fmpz_holonomic_print(op, "n", "Sn"); printf("\n\n");
            printf("start = %lu, n = %lu\n", start, n);

            fmpz_mat_print_pretty(M1); printf("\n\n");
            fmpz_print(Q1); printf("\n\n");

            nmod_mat_print_pretty(M2); printf("\n\n");
            printf("%lu\n\n", Q2);

            abort();
        }

        fmpz_mat_clear(M1);
        fmpz_clear(Q1);
        nmod_mat_clear(M2);
        nmod_mat_clear(L);
        nmod_mat_clear(R);

        fmpz_holonomic_clear(op);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

