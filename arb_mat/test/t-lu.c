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

int fmpq_mat_is_invertible(const fmpq_mat_t A)
{
    int r;
    fmpq_t t;
    fmpq_init(t);
    fmpq_mat_det(t, A);
    r = !fmpq_is_zero(t);
    fmpq_clear(t);
    return r;
}

int main()
{
    long iter;
    flint_rand_t state;

    printf("lu....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        fmpq_mat_t Q;
        arb_mat_t A, LU, P, L, U, T;
        long i, j, n, qbits, prec, *perm;
        int q_invertible, r_invertible;

        n = n_randint(state, 8);
        qbits = 1 + n_randint(state, 100);
        prec = 2 + n_randint(state, 202);

        fmpq_mat_init(Q, n, n);
        arb_mat_init(A, n, n);
        arb_mat_init(LU, n, n);
        arb_mat_init(P, n, n);
        arb_mat_init(L, n, n);
        arb_mat_init(U, n, n);
        arb_mat_init(T, n, n);
        perm = _perm_init(n);

        fmpq_mat_randtest(Q, state, qbits);
        q_invertible = fmpq_mat_is_invertible(Q);

        if (!q_invertible)
        {
            arb_mat_set_fmpq_mat(A, Q, prec);
            r_invertible = arb_mat_lu(perm, LU, A, prec);
            if (r_invertible)
            {
                printf("FAIL: matrix is singular over Q but not over R\n");
                printf("n = %ld, prec = %ld\n", n, prec);
                printf("\n");

                printf("Q = \n"); fmpq_mat_print(Q); printf("\n\n");
                printf("A = \n"); arb_mat_printd(A, 15); printf("\n\n");
                printf("LU = \n"); arb_mat_printd(LU, 15); printf("\n\n");
            }
        }
        else
        {
            /* now this must converge */
            while (1)
            {
                arb_mat_set_fmpq_mat(A, Q, prec);
                r_invertible = arb_mat_lu(perm, LU, A, prec);
                if (r_invertible)
                {
                    break;
                }
                else
                {
                    if (prec > 10000)
                    {
                        printf("FAIL: failed to converge at 10000 bits\n");
                        abort();
                    }
                    prec *= 2;
                }
            }

            arb_mat_one(L);
            for (i = 0; i < n; i++)
                for (j = 0; j < i; j++)
                    arb_set(arb_mat_entry(L, i, j),
                        arb_mat_entry(LU, i, j));

            for (i = 0; i < n; i++)
                for (j = i; j < n; j++)
                    arb_set(arb_mat_entry(U, i, j),
                        arb_mat_entry(LU, i, j));

            for (i = 0; i < n; i++)
                arb_one(arb_mat_entry(P, perm[i], i));

            arb_mat_mul(T, P, L, prec);
            arb_mat_mul(T, T, U, prec);

            if (!arb_mat_contains_fmpq_mat(T, Q))
            {
                printf("FAIL (containment, iter = %ld)\n", iter);
                printf("n = %ld, prec = %ld\n", n, prec);
                printf("\n");

                printf("Q = \n"); fmpq_mat_print(Q); printf("\n\n");
                printf("A = \n"); arb_mat_printd(A, 15); printf("\n\n");
                printf("LU = \n"); arb_mat_printd(LU, 15); printf("\n\n");
                printf("L = \n"); arb_mat_printd(L, 15); printf("\n\n");
                printf("U = \n"); arb_mat_printd(U, 15); printf("\n\n");
                printf("P*L*U = \n"); arb_mat_printd(T, 15); printf("\n\n");

                abort();
            }
        }

        fmpq_mat_clear(Q);
        arb_mat_clear(A);
        arb_mat_clear(LU);
        arb_mat_clear(P);
        arb_mat_clear(L);
        arb_mat_clear(U);
        arb_mat_clear(T);
        _perm_clear(perm);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
