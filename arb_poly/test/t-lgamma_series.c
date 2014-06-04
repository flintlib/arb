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

    Copyright (C) 2012, 2013 Fredrik Johansson

******************************************************************************/

#include "arb_poly.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("lgamma_series....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        long m, n1, n2, qbits, rbits1, rbits2, rbits3;
        fmpq_poly_t A;
        arb_poly_t a, b, c, d;

        qbits = 2 + n_randint(state, 200);
        rbits1 = 2 + n_randint(state, 400);
        rbits2 = 2 + n_randint(state, 400);
        rbits3 = 2 + n_randint(state, 400);

        m = 1 + n_randint(state, 30);
        n1 = 1 + n_randint(state, 30);
        n2 = 1 + n_randint(state, 30);

        fmpq_poly_init(A);
        arb_poly_init(a);
        arb_poly_init(b);
        arb_poly_init(c);
        arb_poly_init(d);

        fmpq_poly_randtest_not_zero(A, state, m, qbits);
        fmpz_abs(A->coeffs, A->coeffs);
        arb_poly_set_fmpq_poly(a, A, rbits1);

        arb_poly_lgamma_series(b, a, n1, rbits2);
        arb_poly_lgamma_series(c, a, n2, rbits3);

        arb_poly_set(d, b);
        arb_poly_truncate(d, FLINT_MIN(n1, n2));
        arb_poly_truncate(c, FLINT_MIN(n1, n2));

        if (!arb_poly_overlaps(c, d))
        {
            printf("FAIL\n\n");
            printf("n1 = %ld, n2 = %ld, bits2 = %ld, bits3 = %ld\n", n1, n2, rbits2, rbits3);

            printf("A = "); fmpq_poly_print(A); printf("\n\n");
            printf("a = "); arb_poly_printd(a, 15); printf("\n\n");
            printf("b = "); arb_poly_printd(b, 15); printf("\n\n");
            printf("c = "); arb_poly_printd(c, 15); printf("\n\n");

            abort();
        }

        /* check loggamma(a) + log(a) = loggamma(a+1) */
        arb_poly_log_series(c, a, n1, rbits2);
        arb_poly_add(c, b, c, rbits2);

        arb_poly_set(d, a);
        arb_add_ui(d->coeffs, d->coeffs, 1, rbits2);
        arb_poly_lgamma_series(d, d, n1, rbits2);

        if (!arb_poly_overlaps(c, d))
        {
            printf("FAIL (functional equation)\n\n");

            printf("A = "); fmpq_poly_print(A); printf("\n\n");
            printf("a = "); arb_poly_printd(a, 15); printf("\n\n");
            printf("b = "); arb_poly_printd(b, 15); printf("\n\n");
            printf("c = "); arb_poly_printd(c, 15); printf("\n\n");
            printf("d = "); arb_poly_printd(d, 15); printf("\n\n");

            abort();
        }

        arb_poly_lgamma_series(a, a, n1, rbits2);
        if (!arb_poly_overlaps(a, b))
        {
            printf("FAIL (aliasing)\n\n");
            abort();
        }

        fmpq_poly_clear(A);
        arb_poly_clear(a);
        arb_poly_clear(b);
        arb_poly_clear(c);
        arb_poly_clear(d);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

