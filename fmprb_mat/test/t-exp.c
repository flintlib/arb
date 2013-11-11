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

#include "fmprb_mat.h"

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
        fmprb_mat_t A, E, F, EF, G;
        fmpq_mat_t Q;
        fmprb_t c, d;
        long n, qbits, prec;

        n = n_randint(state, 5);
        qbits = 2 + n_randint(state, 300);
        prec = 2 + n_randint(state, 300);

        fmprb_init(c);
        fmprb_init(d);
        fmpq_mat_init(Q, n, n);
        fmprb_mat_init(A, n, n);
        fmprb_mat_init(E, n, n);
        fmprb_mat_init(F, n, n);
        fmprb_mat_init(EF, n, n);
        fmprb_mat_init(G, n, n);

        fmpq_mat_randtest(Q, state, qbits);
        fmprb_mat_set_fmpq_mat(A, Q, prec);

        fmprb_mat_exp(E, A, prec);

        fmprb_randtest(c, state, prec, 10);
        fmprb_mat_scalar_mul_fmprb(F, A, c, prec);
        fmprb_mat_exp(F, F, prec);

        fmprb_add_ui(d, c, 1, prec);
        fmprb_mat_scalar_mul_fmprb(G, A, d, prec);
        fmprb_mat_exp(G, G, prec);

        fmprb_mat_mul(EF, E, F, prec);

        if (!fmprb_mat_overlaps(EF, G))
        {
            printf("FAIL\n\n");
            printf("n = %ld, prec = %ld\n", n, prec);

            printf("c = \n"); fmprb_printd(c, 15); printf("\n\n");

            printf("A = \n"); fmprb_mat_printd(A, 15); printf("\n\n");
            printf("E   = \n"); fmprb_mat_printd(E, 15); printf("\n\n");
            printf("F   = \n"); fmprb_mat_printd(F, 15); printf("\n\n");
            printf("E*F = \n"); fmprb_mat_printd(EF, 15); printf("\n\n");
            printf("G   = \n"); fmprb_mat_printd(G, 15); printf("\n\n");

            abort();
        }

        fmprb_clear(c);
        fmprb_clear(d);
        fmpq_mat_clear(Q);
        fmprb_mat_clear(A);
        fmprb_mat_clear(E);
        fmprb_mat_clear(F);
        fmprb_mat_clear(EF);
        fmprb_mat_clear(G);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

