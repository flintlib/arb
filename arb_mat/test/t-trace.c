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

    Copyright (C) 2015 Arb authors

******************************************************************************/

#include "arb_mat.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("trace....");
    fflush(stdout);

    flint_randinit(state);

    /* check that the arb trace contains the fmpq trace */
    for (iter = 0; iter < 10000; iter++)
    {
        fmpq_mat_t Q;
        fmpq_t Qtrace;
        arb_mat_t A;
        arb_t Atrace;
        slong n, qbits, prec;

        n = n_randint(state, 8);
        qbits = 1 + n_randint(state, 100);
        prec = 2 + n_randint(state, 200);

        fmpq_mat_init(Q, n, n);
        fmpq_init(Qtrace);

        arb_mat_init(A, n, n);
        arb_init(Atrace);

        fmpq_mat_randtest(Q, state, qbits);
        fmpq_mat_trace(Qtrace, Q);

        arb_mat_set_fmpq_mat(A, Q, prec);
        arb_mat_trace(Atrace, A, prec);

        if (!arb_contains_fmpq(Atrace, Qtrace))
        {
            flint_printf("FAIL (containment, iter = %wd)\n", iter);
            flint_printf("n = %wd, prec = %wd\n", n, prec);
            flint_printf("\n");

            flint_printf("Q = \n"); fmpq_mat_print(Q); flint_printf("\n\n");
            flint_printf("Qtrace = \n"); fmpq_print(Qtrace); flint_printf("\n\n");

            flint_printf("A = \n"); arb_mat_printd(A, 15); flint_printf("\n\n");
            flint_printf("Atrace = \n"); arb_printd(Atrace, 15); flint_printf("\n\n");
            flint_printf("Atrace = \n"); arb_print(Atrace); flint_printf("\n\n");

            abort();
        }

        fmpq_mat_clear(Q);
        fmpq_clear(Qtrace);
        arb_mat_clear(A);
        arb_clear(Atrace);
    }

    /* check trace(A*B) = trace(B*A) */
    for (iter = 0; iter < 10000; iter++)
    {
        slong m, n, prec;
        arb_mat_t a, b, ab, ba;
        arb_t trab, trba;

        prec = 2 + n_randint(state, 200);

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        arb_mat_init(a, m, n);
        arb_mat_init(b, n, m);
        arb_mat_init(ab, m, m);
        arb_mat_init(ba, n, n);

        arb_init(trab);
        arb_init(trba);

        arb_mat_randtest(a, state, 2 + n_randint(state, 100), 10);
        arb_mat_randtest(b, state, 2 + n_randint(state, 100), 10);

        arb_mat_mul(ab, a, b, prec);
        arb_mat_mul(ba, b, a, prec);

        arb_mat_trace(trab, ab, prec);
        arb_mat_trace(trba, ba, prec);

        if (!arb_overlaps(trab, trba))
        {
            flint_printf("FAIL (overlap, iter = %wd)\n", iter);
            flint_printf("m = %wd, n = %wd, prec = %wd\n", m, n, prec);
            flint_printf("\n");

            flint_printf("a = \n"); arb_mat_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = \n"); arb_mat_printd(b, 15); flint_printf("\n\n");
            flint_printf("ab = \n"); arb_mat_printd(ab, 15); flint_printf("\n\n");
            flint_printf("ba = \n"); arb_mat_printd(ba, 15); flint_printf("\n\n");

            flint_printf("trace(ab) = \n"); arb_printd(trab, 15); flint_printf("\n\n");
            flint_printf("trace(ba) = \n"); arb_printd(trba, 15); flint_printf("\n\n");
            abort();
        }

        arb_clear(trab);
        arb_clear(trba);

        arb_mat_clear(a);
        arb_mat_clear(b);
        arb_mat_clear(ab);
        arb_mat_clear(ba);
    }

    /* check trace(A^T A) = frobenius_norm(A)^2 */
    for (iter = 0; iter < 10000; iter++)
    {
        slong m, n, prec;
        arb_mat_t A, AT, ATA;
        arb_t t;
        mag_t low, frobenius;

        prec = 2 + n_randint(state, 200);

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        arb_mat_init(A, m, n);
        arb_mat_randtest(A, state, 2 + n_randint(state, 100), 10);

        arb_mat_init(AT, n, m);
        arb_mat_transpose(AT, A);

        arb_mat_init(ATA, n, n);
        arb_mat_mul(ATA, AT, A, prec);

        arb_init(t);
        arb_mat_trace(t, ATA, prec);
        arb_sqrt(t, t, prec);

        mag_init(low);
        arb_get_mag_lower(low, t);

        mag_init(frobenius);
        arb_mat_bound_frobenius_norm(frobenius, A);

        if (mag_cmp(low, frobenius) > 0)
        {
            flint_printf("FAIL (frobenius norm)\n", iter);
            flint_printf("m = %wd, n = %wd, prec = %wd\n", m, n, prec);
            flint_printf("\n");

            flint_printf("A = \n"); arb_mat_printd(A, 15); flint_printf("\n\n");

            flint_printf("lower(sqrt(trace(A^T A))) = \n");
            mag_printd(low, 15); flint_printf("\n\n");

            flint_printf("upper(frobenius_norm(A)) = \n");
            mag_printd(frobenius, 15); flint_printf("\n\n");

            abort();
        }

        arb_clear(t);

        mag_clear(low);
        mag_clear(frobenius);

        arb_mat_clear(A);
        arb_mat_clear(AT);
        arb_mat_clear(ATA);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
