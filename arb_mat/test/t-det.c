/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("det....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000 * arb_test_multiplier(); iter++)
    {
        fmpq_mat_t Q;
        fmpq_t Qdet;
        arb_mat_t A;
        arb_t Adet;
        slong n, qbits, prec;

        n = n_randint(state, 8);
        qbits = 1 + n_randint(state, 100);
        prec = 2 + n_randint(state, 200);

        fmpq_mat_init(Q, n, n);
        fmpq_init(Qdet);

        arb_mat_init(A, n, n);
        arb_init(Adet);

        fmpq_mat_randtest(Q, state, qbits);
        fmpq_mat_det(Qdet, Q);

        arb_mat_set_fmpq_mat(A, Q, prec);
        arb_mat_det(Adet, A, prec);

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

            flint_abort();
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
