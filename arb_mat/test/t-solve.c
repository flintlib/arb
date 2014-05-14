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

#include "arb_mat.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("solve....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        fmpq_mat_t Q, QX, QB;
        arb_mat_t A, X, B;
        long n, m, qbits, prec;
        int q_invertible, r_invertible, r_invertible2;

        n = n_randint(state, 8);
        m = n_randint(state, 8);
        qbits = 1 + n_randint(state, 30);
        prec = 2 + n_randint(state, 200);

        fmpq_mat_init(Q, n, n);
        fmpq_mat_init(QX, n, m);
        fmpq_mat_init(QB, n, m);

        arb_mat_init(A, n, n);
        arb_mat_init(X, n, m);
        arb_mat_init(B, n, m);

        fmpq_mat_randtest(Q, state, qbits);
        fmpq_mat_randtest(QB, state, qbits);

        q_invertible = fmpq_mat_solve_fraction_free(QX, Q, QB);

        if (!q_invertible)
        {
            arb_mat_set_fmpq_mat(A, Q, prec);
            r_invertible = arb_mat_solve(X, A, B, prec);
            if (r_invertible)
            {
                printf("FAIL: matrix is singular over Q but not over R\n");
                printf("n = %ld, prec = %ld\n", n, prec);
                printf("\n");

                printf("Q = \n"); fmpq_mat_print(Q); printf("\n\n");
                printf("QX = \n"); fmpq_mat_print(QX); printf("\n\n");
                printf("QB = \n"); fmpq_mat_print(QB); printf("\n\n");
                printf("A = \n"); arb_mat_printd(A, 15); printf("\n\n");
                abort();
            }
        }
        else
        {
            /* now this must converge */
            while (1)
            {
                arb_mat_set_fmpq_mat(A, Q, prec);
                arb_mat_set_fmpq_mat(B, QB, prec);

                r_invertible = arb_mat_solve(X, A, B, prec);
                if (r_invertible)
                {
                    break;
                }
                else
                {
                    if (prec > 10000)
                    {
                        printf("FAIL: failed to converge at 10000 bits\n");
                        printf("Q = \n"); fmpq_mat_print(Q); printf("\n\n");
                        printf("QX = \n"); fmpq_mat_print(QX); printf("\n\n");
                        printf("QB = \n"); fmpq_mat_print(QB); printf("\n\n");
                        printf("A = \n"); arb_mat_printd(A, 15); printf("\n\n");
                        abort();
                    }
                    prec *= 2;
                }
            }

            if (!arb_mat_contains_fmpq_mat(X, QX))
            {
                printf("FAIL (containment, iter = %ld)\n", iter);
                printf("n = %ld, prec = %ld\n", n, prec);
                printf("\n");

                printf("Q = \n"); fmpq_mat_print(Q); printf("\n\n");
                printf("QB = \n"); fmpq_mat_print(QB); printf("\n\n");
                printf("QX = \n"); fmpq_mat_print(QX); printf("\n\n");

                printf("A = \n"); arb_mat_printd(A, 15); printf("\n\n");
                printf("B = \n"); arb_mat_printd(B, 15); printf("\n\n");
                printf("X = \n"); arb_mat_printd(X, 15); printf("\n\n");

                abort();
            }

            /* test aliasing */
            r_invertible2 = arb_mat_solve(B, A, B, prec);
            if (!arb_mat_equal(X, B) || r_invertible != r_invertible2)
            {
                printf("FAIL (aliasing)\n");
                printf("A = \n"); arb_mat_printd(A, 15); printf("\n\n");
                printf("B = \n"); arb_mat_printd(B, 15); printf("\n\n");
                printf("X = \n"); arb_mat_printd(X, 15); printf("\n\n");
                abort();
            }
        }

        fmpq_mat_clear(Q);
        fmpq_mat_clear(QB);
        fmpq_mat_clear(QX);
        arb_mat_clear(A);
        arb_mat_clear(B);
        arb_mat_clear(X);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
