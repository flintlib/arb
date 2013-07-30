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

    printf("forward_fmprb_mat....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        fmpz_holonomic_t op;
        fmpz_mat_t M1;
        fmpz_t Q1;
        fmprb_mat_t M2;
        fmprb_t Q2;
        long start, n, r, prec;

        fmpz_holonomic_init(op);

        fmpz_holonomic_randtest(op, state, 5, 5, 10);

        r = fmpz_holonomic_order(op);
        start = n_randint(state, 20);
        n = n_randint(state, 20);
        prec = 2 + n_randint(state, 200);

        fmpz_mat_init(M1, r, r);
        fmpz_init(Q1);
        fmprb_mat_init(M2, r, r);
        fmprb_init(Q2);

        fmpz_holonomic_forward_fmpz_mat(M1, Q1, op, start, n);
        fmpz_holonomic_forward_fmprb_mat(M2, Q2, op, start, n, prec);

        if (!fmprb_mat_contains_fmpz_mat(M2, M1) ||
            !fmprb_contains_fmpz(Q2, Q1))
        {
            printf("FAIL (containment)\n");

            fmpz_holonomic_print(op, "n", "Sn"); printf("\n\n");

            fmpz_mat_print_pretty(M1); printf("\n\n");
            fmpz_print(Q1); printf("\n\n");

            fmprb_mat_printd(M2, 10); printf("\n\n");
            fmprb_printd(Q2, 10); printf("\n\n");

            abort();
        }

        fmpz_mat_clear(M1);
        fmpz_clear(Q1);
        fmprb_mat_clear(M2);
        fmprb_clear(Q2);

        fmpz_holonomic_clear(op);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

