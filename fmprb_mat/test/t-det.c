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

#include "fmprb_mat.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("det....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        fmpq_mat_t Q;
        fmpq_t Qdet;
        fmprb_mat_t A;
        fmprb_t Adet;
        long n, qbits, prec;

        n = n_randint(state, 8);
        qbits = 1 + n_randint(state, 100);
        prec = 2 + n_randint(state, 200);

        fmpq_mat_init(Q, n, n);
        fmpq_init(Qdet);

        fmprb_mat_init(A, n, n);
        fmprb_init(Adet);

        fmpq_mat_randtest(Q, state, qbits);
        fmpq_mat_det(Qdet, Q);

        fmprb_mat_set_fmpq_mat(A, Q, prec);
        fmprb_mat_det(Adet, A, prec);

        if (!fmprb_contains_fmpq(Adet, Qdet))
        {
            printf("FAIL (containment, iter = %ld)\n", iter);
            printf("n = %ld, prec = %ld\n", n, prec);
            printf("\n");

            printf("Q = \n"); fmpq_mat_print(Q); printf("\n\n");
            printf("Qdet = \n"); fmpq_print(Qdet); printf("\n\n");

            printf("A = \n"); fmprb_mat_printd(A, 15); printf("\n\n");
            printf("Adet = \n"); fmprb_printd(Adet, 15); printf("\n\n");
            printf("Adet = \n"); fmprb_print(Adet); printf("\n\n");

            abort();
        }

        fmpq_mat_clear(Q);
        fmpq_clear(Qdet);
        fmprb_mat_clear(A);
        fmprb_clear(Adet);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
