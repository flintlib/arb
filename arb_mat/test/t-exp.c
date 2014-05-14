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

#include "arb_mat.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("exp....");
    fflush(stdout);

    flint_randinit(state);

    /* check exp(A)*exp(c*A) = exp((1+c)*A) */
    for (iter = 0; iter < 1000; iter++)
    {
        arb_mat_t A, E, F, EF, G;
        fmpq_mat_t Q;
        arb_t c, d;
        long n, qbits, prec;

        n = n_randint(state, 5);
        qbits = 2 + n_randint(state, 300);
        prec = 2 + n_randint(state, 300);

        arb_init(c);
        arb_init(d);
        fmpq_mat_init(Q, n, n);
        arb_mat_init(A, n, n);
        arb_mat_init(E, n, n);
        arb_mat_init(F, n, n);
        arb_mat_init(EF, n, n);
        arb_mat_init(G, n, n);

        fmpq_mat_randtest(Q, state, qbits);
        arb_mat_set_fmpq_mat(A, Q, prec);

        arb_mat_exp(E, A, prec);

        arb_randtest(c, state, prec, 10);
        arb_mat_scalar_mul_arb(F, A, c, prec);
        arb_mat_exp(F, F, prec);

        arb_add_ui(d, c, 1, prec);
        arb_mat_scalar_mul_arb(G, A, d, prec);
        arb_mat_exp(G, G, prec);

        arb_mat_mul(EF, E, F, prec);

        if (!arb_mat_overlaps(EF, G))
        {
            printf("FAIL\n\n");
            printf("n = %ld, prec = %ld\n", n, prec);

            printf("c = \n"); arb_printd(c, 15); printf("\n\n");

            printf("A = \n"); arb_mat_printd(A, 15); printf("\n\n");
            printf("E   = \n"); arb_mat_printd(E, 15); printf("\n\n");
            printf("F   = \n"); arb_mat_printd(F, 15); printf("\n\n");
            printf("E*F = \n"); arb_mat_printd(EF, 15); printf("\n\n");
            printf("G   = \n"); arb_mat_printd(G, 15); printf("\n\n");

            abort();
        }

        arb_clear(c);
        arb_clear(d);
        fmpq_mat_clear(Q);
        arb_mat_clear(A);
        arb_mat_clear(E);
        arb_mat_clear(F);
        arb_mat_clear(EF);
        arb_mat_clear(G);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

