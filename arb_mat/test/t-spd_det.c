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

    Copyright (C) 2016 Arb authors

******************************************************************************/

#include "arb_mat.h"

void
_fmpq_mat_randtest_positive_semidefinite(fmpq_mat_t mat, flint_rand_t state, mp_bitcnt_t bits)
{
    slong n;
    fmpq_mat_t R, RT;
    if (!fmpq_mat_is_square(mat)) abort(); /* assert */
    n = fmpq_mat_nrows(mat);
    fmpq_mat_init(R, n, n);
    fmpq_mat_init(RT, n, n);
    fmpq_mat_randtest(R, state, bits);
    fmpq_mat_transpose(RT, R);
    fmpq_mat_mul(mat, R, RT);
    fmpq_mat_clear(R);
    fmpq_mat_clear(RT);
}

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("spd_det....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        fmpq_mat_t Q;
        fmpq_t Qdet;
        arb_mat_t A;
        arb_t Adet;
        slong n, qbits, prec;
        int r_invertible, q_invertible;

        n = n_randint(state, 8);
        qbits = 1 + n_randint(state, 100);
        prec = 2 + n_randint(state, 200);

        fmpq_mat_init(Q, n, n);
        fmpq_init(Qdet);

        arb_mat_init(A, n, n);
        arb_init(Adet);

        _fmpq_mat_randtest_positive_semidefinite(Q, state, qbits);
        fmpq_mat_det(Qdet, Q);
        q_invertible = !fmpq_is_zero(Qdet);

        if (!q_invertible)
        {
            arb_mat_set_fmpq_mat(A, Q, prec);
            r_invertible = arb_mat_spd_det(Adet, A, prec);
            if (r_invertible)
            {
                flint_printf("FAIL: matrix is singular over Q but not over R\n");
                flint_printf("n = %wd, prec = %wd\n", n, prec);
                flint_printf("\n");

                flint_printf("Q = \n"); fmpq_mat_print(Q); flint_printf("\n\n");
                flint_printf("A = \n"); arb_mat_printd(A, 15); flint_printf("\n\n");
                abort();
            }
        }
        else
        {
            /* now this must converge */
            while (1)
            {
                arb_mat_set_fmpq_mat(A, Q, prec);
                r_invertible = arb_mat_spd_det(Adet, A, prec);

                if (r_invertible)
                {
                    break;
                }
                else
                {
                    if (prec > 10000)
                    {
                        flint_printf("FAIL: failed to converge at 10000 bits\n");
                        flint_printf("Q = \n"); fmpq_mat_print(Q); flint_printf("\n\n");
                        flint_printf("A = \n"); arb_mat_printd(A, 15); flint_printf("\n\n");
                        abort();
                    }
                    prec *= 2;
                }
            }

            if (!arb_contains_fmpq(Adet, Qdet))
            {
                flint_printf("FAIL (containment, iter = %wd)\n", iter);
                flint_printf("n = %wd, prec = %wd\n", n, prec);
                flint_printf("\n");

                flint_printf("Q = \n"); fmpq_mat_print(Q); flint_printf("\n\n");
                flint_printf("Qdet = \n"); fmpq_print(Qdet); flint_printf("\n\n");

                flint_printf("A = \n"); arb_mat_printd(A, 15); flint_printf("\n\n");
                flint_printf("Adet = \n"); arb_printd(Adet, 15); flint_printf("\n\n");
                flint_printf("Adet = \n"); arb_print(Adet); flint_printf("\n\n");
                abort();
            }

        }

        fmpq_mat_clear(Q);
        fmpq_clear(Qdet);
        arb_mat_clear(A);
        arb_clear(Adet);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
