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

#include "acb_mat.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("exp....");
    fflush(stdout);

    flint_randinit(state);

    /* check exp(A)*exp(c*A) = exp((1+c)*A) */
    for (iter = 0; iter < 500; iter++)
    {
        acb_mat_t A, E, F, EF, G;
        fmpq_mat_t Q;
        acb_t c, d;
        long n, qbits, prec;

        n = n_randint(state, 5);
        qbits = 2 + n_randint(state, 300);
        prec = 2 + n_randint(state, 300);

        acb_init(c);
        acb_init(d);
        fmpq_mat_init(Q, n, n);
        acb_mat_init(A, n, n);
        acb_mat_init(E, n, n);
        acb_mat_init(F, n, n);
        acb_mat_init(EF, n, n);
        acb_mat_init(G, n, n);

        fmpq_mat_randtest(Q, state, qbits);
        acb_mat_set_fmpq_mat(A, Q, prec);

        acb_mat_exp(E, A, prec);

        acb_randtest(c, state, prec, 10);
        acb_mat_scalar_mul_acb(F, A, c, prec);
        acb_mat_exp(F, F, prec);

        acb_add_ui(d, c, 1, prec);
        acb_mat_scalar_mul_acb(G, A, d, prec);
        acb_mat_exp(G, G, prec);

        acb_mat_mul(EF, E, F, prec);

        if (!acb_mat_overlaps(EF, G))
        {
            printf("FAIL\n\n");
            printf("n = %ld, prec = %ld\n", n, prec);

            printf("c = \n"); acb_printd(c, 15); printf("\n\n");

            printf("A = \n"); acb_mat_printd(A, 15); printf("\n\n");
            printf("E   = \n"); acb_mat_printd(E, 15); printf("\n\n");
            printf("F   = \n"); acb_mat_printd(F, 15); printf("\n\n");
            printf("E*F = \n"); acb_mat_printd(EF, 15); printf("\n\n");
            printf("G   = \n"); acb_mat_printd(G, 15); printf("\n\n");

            abort();
        }

        acb_clear(c);
        acb_clear(d);
        fmpq_mat_clear(Q);
        acb_mat_clear(A);
        acb_mat_clear(E);
        acb_mat_clear(F);
        acb_mat_clear(EF);
        acb_mat_clear(G);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

