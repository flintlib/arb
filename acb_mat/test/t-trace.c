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

#include "acb_mat.h"

static void
_acb_mat_conjugate_transpose(acb_mat_t B, const acb_mat_t A)
{
    slong i, j;
    acb_mat_transpose(B, A);
    for (i = 0; i < acb_mat_nrows(B); i++)
        for (j = 0; j < acb_mat_ncols(B); j++)
            acb_conj(acb_mat_entry(B, i, j), acb_mat_entry(B, i, j));
}

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("trace....");
    fflush(stdout);

    flint_randinit(state);

    /* check that the acb trace contains the fmpq trace */
    for (iter = 0; iter < 10000; iter++)
    {
        fmpq_mat_t Q;
        fmpq_t Qtrace;
        acb_mat_t A;
        acb_t Atrace;
        slong n, qbits, prec;

        n = n_randint(state, 8);
        qbits = 1 + n_randint(state, 100);
        prec = 2 + n_randint(state, 200);

        fmpq_mat_init(Q, n, n);
        fmpq_init(Qtrace);

        acb_mat_init(A, n, n);
        acb_init(Atrace);

        fmpq_mat_randtest(Q, state, qbits);
        fmpq_mat_trace(Qtrace, Q);

        acb_mat_set_fmpq_mat(A, Q, prec);
        acb_mat_trace(Atrace, A, prec);

        if (!acb_contains_fmpq(Atrace, Qtrace))
        {
            flint_printf("FAIL (containment, iter = %wd)\n", iter);
            flint_printf("n = %wd, prec = %wd\n", n, prec);
            flint_printf("\n");

            flint_printf("Q = \n"); fmpq_mat_print(Q); flint_printf("\n\n");
            flint_printf("Qtrace = \n"); fmpq_print(Qtrace); flint_printf("\n\n");

            flint_printf("A = \n"); acb_mat_printd(A, 15); flint_printf("\n\n");
            flint_printf("Atrace = \n"); acb_printd(Atrace, 15); flint_printf("\n\n");
            flint_printf("Atrace = \n"); acb_print(Atrace); flint_printf("\n\n");

            abort();
        }

        fmpq_mat_clear(Q);
        fmpq_clear(Qtrace);
        acb_mat_clear(A);
        acb_clear(Atrace);
    }

    /* check trace(A*B) = trace(B*A) */
    for (iter = 0; iter < 10000; iter++)
    {
        slong m, n, prec;
        acb_mat_t a, b, ab, ba;
        acb_t trab, trba;

        prec = 2 + n_randint(state, 200);

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        acb_mat_init(a, m, n);
        acb_mat_init(b, n, m);
        acb_mat_init(ab, m, m);
        acb_mat_init(ba, n, n);

        acb_init(trab);
        acb_init(trba);

        acb_mat_randtest(a, state, 2 + n_randint(state, 100), 10);
        acb_mat_randtest(b, state, 2 + n_randint(state, 100), 10);

        acb_mat_mul(ab, a, b, prec);
        acb_mat_mul(ba, b, a, prec);

        acb_mat_trace(trab, ab, prec);
        acb_mat_trace(trba, ba, prec);

        if (!acb_overlaps(trab, trba))
        {
            flint_printf("FAIL (overlap, iter = %wd)\n", iter);
            flint_printf("m = %wd, n = %wd, prec = %wd\n", m, n, prec);
            flint_printf("\n");

            flint_printf("a = \n"); acb_mat_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = \n"); acb_mat_printd(b, 15); flint_printf("\n\n");
            flint_printf("ab = \n"); acb_mat_printd(ab, 15); flint_printf("\n\n");
            flint_printf("ba = \n"); acb_mat_printd(ba, 15); flint_printf("\n\n");

            flint_printf("trace(ab) = \n"); acb_printd(trab, 15); flint_printf("\n\n");
            flint_printf("trace(ba) = \n"); acb_printd(trba, 15); flint_printf("\n\n");
        }

        acb_clear(trab);
        acb_clear(trba);

        acb_mat_clear(a);
        acb_mat_clear(b);
        acb_mat_clear(ab);
        acb_mat_clear(ba);
    }

    /* check trace(A^H A) = frobenius_norm(A)^2 */
    for (iter = 0; iter < 10000; iter++)
    {
        slong m, n, prec;
        acb_mat_t A, AH, AHA;
        acb_t t;
        mag_t low, fro;

        prec = 2 + n_randint(state, 200);

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        acb_mat_init(A, m, n);
        acb_mat_randtest(A, state, 2 + n_randint(state, 100), 10);

        acb_mat_init(AH, n, m);
        _acb_mat_conjugate_transpose(AH, A);

        acb_mat_init(AHA, n, n);
        acb_mat_mul(AHA, AH, A, prec);

        acb_init(t);
        acb_mat_trace(t, AHA, prec);
        acb_sqrt(t, t, prec);

        mag_init(low);
        acb_get_mag_lower(low, t);

        mag_init(fro);
        acb_mat_bound_fro_norm(fro, A);

        if (mag_cmp(low, fro) > 0)
        {
            flint_printf("FAIL (frobenius norm)\n", iter);
            flint_printf("m = %wd, n = %wd, prec = %wd\n", m, n, prec);
            flint_printf("\n");

            flint_printf("A = \n"); acb_mat_printd(A, 15); flint_printf("\n\n");

            flint_printf("lower(sqrt(trace(A^H A))) = \n");
            mag_printd(low, 15); flint_printf("\n\n");

            flint_printf("upper(frobenius_norm(A)) = \n");
            mag_printd(fro, 15); flint_printf("\n\n");

            abort();
        }

        acb_clear(t);

        mag_clear(low);
        mag_clear(fro);

        acb_mat_clear(A);
        acb_mat_clear(AH);
        acb_mat_clear(AHA);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
