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

#include "fmpcb_mat.h"

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
        fmpcb_mat_t A, E, F, EF, G;
        fmpq_mat_t Q;
        fmpcb_t c, d;
        long n, qbits, prec;

        n = n_randint(state, 5);
        qbits = 2 + n_randint(state, 300);
        prec = 2 + n_randint(state, 300);

        fmpcb_init(c);
        fmpcb_init(d);
        fmpq_mat_init(Q, n, n);
        fmpcb_mat_init(A, n, n);
        fmpcb_mat_init(E, n, n);
        fmpcb_mat_init(F, n, n);
        fmpcb_mat_init(EF, n, n);
        fmpcb_mat_init(G, n, n);

        fmpq_mat_randtest(Q, state, qbits);
        fmpcb_mat_set_fmpq_mat(A, Q, prec);

        fmpcb_mat_exp(E, A, prec);

        fmpcb_randtest(c, state, prec, 10);
        fmpcb_mat_scalar_mul_fmpcb(F, A, c, prec);
        fmpcb_mat_exp(F, F, prec);

        fmpcb_add_ui(d, c, 1, prec);
        fmpcb_mat_scalar_mul_fmpcb(G, A, d, prec);
        fmpcb_mat_exp(G, G, prec);

        fmpcb_mat_mul(EF, E, F, prec);

        if (!fmpcb_mat_overlaps(EF, G))
        {
            printf("FAIL\n\n");
            printf("n = %ld, prec = %ld\n", n, prec);

            printf("c = \n"); fmpcb_printd(c, 15); printf("\n\n");

            printf("A = \n"); fmpcb_mat_printd(A, 15); printf("\n\n");
            printf("E   = \n"); fmpcb_mat_printd(E, 15); printf("\n\n");
            printf("F   = \n"); fmpcb_mat_printd(F, 15); printf("\n\n");
            printf("E*F = \n"); fmpcb_mat_printd(EF, 15); printf("\n\n");
            printf("G   = \n"); fmpcb_mat_printd(G, 15); printf("\n\n");

            abort();
        }

        fmpcb_clear(c);
        fmpcb_clear(d);
        fmpq_mat_clear(Q);
        fmpcb_mat_clear(A);
        fmpcb_mat_clear(E);
        fmpcb_mat_clear(F);
        fmpcb_mat_clear(EF);
        fmpcb_mat_clear(G);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

