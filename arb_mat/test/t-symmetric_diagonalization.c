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

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("symmetric_diagonalization....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        /* flint_printf("symmetric diagonalization iter %wd\n", iter); */

        slong i, j, n;
        slong prec;
        arb_mat_t D, P, A;

        n = n_randint(state, 8);
        prec = 2 + n_randint(state, 202);

        arb_mat_init(D, n, 1);
        arb_mat_init(P, n, n);
        arb_mat_init(A, n, n);

        arb_mat_randtest(A, state, 2 + n_randint(state, 100), 10);
        for (i = 0; i < n; i++)
            for (j = 0; j < i; j++)
                arb_set(arb_mat_entry(A, i, j), arb_mat_entry(A, j, i));

        /*
        prec = 114;
        arb_one(arb_mat_entry(A, 0, 0));
        arb_zero(arb_mat_entry(A, 1, 1));
        arb_set_d(arb_mat_entry(A, 0, 1), 1e39);
        arb_set_d(arb_mat_entry(A, 1, 0), 1e39);
        */

        /*
        prec = 64;
        arb_mat_zero(A);
        arb_one(arb_mat_entry(A, 0, 0));
        arb_zero(arb_mat_entry(A, 1, 1));
        arf_set_si_2exp_si(arb_midref(arb_mat_entry(A, 0, 1)), 1, 60);
        arf_set_si_2exp_si(arb_midref(arb_mat_entry(A, 1, 0)), 1, 60);
        */

        /*
        prec = 114
            A = 
            [1.29185485839844 +/- 9.8608e-32, 3.34374580571462e+38 +/- 3.3554e+07]
            [3.34374580571462e+38 +/- 3.3554e+07, 1.66605735186693e-13 +/- 0]

        */

        /*
        flint_printf("before diagonalization:\n");
        flint_printf("prec = %wd\n", prec);
        flint_printf("A = \n"); arb_mat_printd(A, 15); flint_printf("\n");
        */

        arb_mat_symmetric_diagonalization(D, P, A, prec);

        /* charpoly should contain zero at eigenvalues of A */
        {
            slong i;
            arb_poly_t f;
            arb_t y;

            arb_poly_init(f);
            arb_init(y);

            arb_mat_charpoly(f, A, prec);

            for (i = 0; i < n; i++)
            {
                arb_srcptr x = arb_mat_entry(D, i, 0);

                /*
                flint_printf("eigenvalue %wd : ", i);
                arb_printd(x, 15); flint_printf("\n");
                */

                arb_poly_evaluate(y, f, x, prec);
                if (!arb_contains_zero(y))
                {
                    flint_printf("FAIL (charpoly)\n");
                    flint_printf("A = \n"); arb_mat_printd(A, 15); flint_printf("\n");
                    flint_printf("D = \n"); arb_mat_printd(D, 15); flint_printf("\n");
                    flint_printf("P = \n"); arb_mat_printd(P, 15); flint_printf("\n");
                    abort();
                }
            }

            arb_poly_clear(f);
            arb_clear(y);
        }

        /* Multiplying out the decomposition should give the original matrix */
        {
            /* A = P * D * P^T */
            slong i, j, k;
            arb_t t;
            arb_mat_t Y;

            arb_init(t);
            arb_mat_init(Y, n, n);

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    for (k = 0; k < n; k++)
                    {
                        arb_mul(t, arb_mat_entry(P, i, k),
                                   arb_mat_entry(P, j, k), prec);
                        arb_addmul(arb_mat_entry(Y, i, j),
                                   arb_mat_entry(D, k, 0), t, prec);
                    }
                }
            }

            if (!arb_mat_contains(Y, A))
            {
                flint_printf("FAIL (decomposition expansion)\n");
                flint_printf("A = \n"); arb_mat_printd(A, 15);
                flint_printf("D = \n"); arb_mat_printd(D, 15);
                flint_printf("P = \n"); arb_mat_printd(P, 15);
                flint_printf("P * D * P^T = \n"); arb_mat_printd(Y, 15);
                abort();
            }

            arb_clear(t);
            arb_mat_clear(Y);
        }

        arb_mat_clear(D);
        arb_mat_clear(P);
        arb_mat_clear(A);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
