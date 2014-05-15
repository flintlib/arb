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

    printf("det....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        fmpq_mat_t Q;
        fmpq_t Qdet;
        acb_mat_t A;
        acb_t Adet, imagunit;
        long n, qbits, prec;
        int imaginary;

        n = n_randint(state, 8);
        qbits = 1 + n_randint(state, 100);
        prec = 2 + n_randint(state, 200);
        imaginary = n_randint(state, 2);

        fmpq_mat_init(Q, n, n);
        fmpq_init(Qdet);

        acb_mat_init(A, n, n);
        acb_init(Adet);
        acb_init(imagunit);

        fmpq_mat_randtest(Q, state, qbits);
        fmpq_mat_det(Qdet, Q);

        acb_mat_set_fmpq_mat(A, Q, prec);

        if (imaginary)
        {
            acb_onei(imagunit);
            acb_mat_scalar_mul_acb(A, A, imagunit, prec);
        }

        acb_mat_det(Adet, A, prec);

        if (imaginary)
        {
            acb_onei(imagunit);
            acb_inv(imagunit, imagunit, prec);
            acb_pow_ui(imagunit, imagunit, n, prec);
            acb_mul(Adet, Adet, imagunit, prec);
        }

        if (!acb_contains_fmpq(Adet, Qdet))
        {
            printf("FAIL (containment, iter = %ld)\n", iter);
            printf("n = %ld, prec = %ld\n", n, prec);
            printf("\n");

            printf("Q = \n"); fmpq_mat_print(Q); printf("\n\n");
            printf("Qdet = \n"); fmpq_print(Qdet); printf("\n\n");

            printf("A = \n"); acb_mat_printd(A, 15); printf("\n\n");
            printf("Adet = \n"); acb_printd(Adet, 15); printf("\n\n");
            printf("Adet = \n"); acb_print(Adet); printf("\n\n");

            abort();
        }

        fmpq_mat_clear(Q);
        fmpq_clear(Qdet);
        acb_mat_clear(A);
        acb_clear(Adet);
        acb_clear(imagunit);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
