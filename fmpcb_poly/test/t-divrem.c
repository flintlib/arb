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

#include "fmpcb_poly.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("divrem....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        long m, n, qbits1, qbits2, rbits1, rbits2, rbits3;
        fmpq_poly_t A, B, Q, R;
        fmpcb_poly_t a, b, q, r;

        qbits1 = 2 + n_randint(state, 200);
        qbits2 = 2 + n_randint(state, 200);
        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);
        rbits3 = 2 + n_randint(state, 200);

        m = 1 + n_randint(state, 20);
        n = 1 + n_randint(state, 20);

        fmpq_poly_init(A);
        fmpq_poly_init(B);
        fmpq_poly_init(Q);
        fmpq_poly_init(R);

        fmpcb_poly_init(a);
        fmpcb_poly_init(b);
        fmpcb_poly_init(q);
        fmpcb_poly_init(r);

        fmpq_poly_randtest(A, state, m, qbits1);
        fmpq_poly_randtest_not_zero(B, state, n, qbits2);

        fmpq_poly_divrem(Q, R, A, B);

        fmpcb_poly_set_fmpq_poly(a, A, rbits1);
        fmpcb_poly_set_fmpq_poly(b, B, rbits2);

        fmpcb_poly_divrem(q, r, a, b, rbits3);

        if (!fmpcb_poly_contains_fmpq_poly(q, Q) ||
             !fmpcb_poly_contains_fmpq_poly(r, R))
        {
            printf("FAIL\n\n");

            printf("A = "); fmpq_poly_print(A); printf("\n\n");
            printf("B = "); fmpq_poly_print(B); printf("\n\n");
            printf("Q = "); fmpq_poly_print(Q); printf("\n\n");
            printf("R = "); fmpq_poly_print(R); printf("\n\n");

            printf("a = "); fmpcb_poly_printd(a, 15); printf("\n\n");
            printf("b = "); fmpcb_poly_printd(b, 15); printf("\n\n");
            printf("q = "); fmpcb_poly_printd(q, 15); printf("\n\n");
            printf("r = "); fmpcb_poly_printd(r, 15); printf("\n\n");

            abort();
        }

        fmpcb_poly_divrem(a, r, a, b, rbits3);
        if (!fmpcb_poly_equal(a, q))
        {
            printf("FAIL (aliasing q, a)\n\n");
        }
        fmpcb_poly_set_fmpq_poly(a, A, rbits1);

        fmpcb_poly_divrem(b, r, a, b, rbits3);
        if (!fmpcb_poly_equal(b, q))
        {
            printf("FAIL (aliasing q, b)\n\n");
            abort();
        }
        fmpcb_poly_set_fmpq_poly(b, B, rbits2);

        fmpcb_poly_divrem(q, a, a, b, rbits3);
        if (!fmpcb_poly_equal(a, r))
        {
            printf("FAIL (aliasing r, a)\n\n");
            abort();
        }
        fmpcb_poly_set_fmpq_poly(a, A, rbits1);

        fmpcb_poly_divrem(q, b, a, b, rbits3);
        if (!fmpcb_poly_equal(b, r))
        {
            printf("FAIL (aliasing r, b)\n\n");
            abort();
        }

        fmpq_poly_clear(A);
        fmpq_poly_clear(B);
        fmpq_poly_clear(Q);
        fmpq_poly_clear(R);

        fmpcb_poly_clear(a);
        fmpcb_poly_clear(b);
        fmpcb_poly_clear(q);
        fmpcb_poly_clear(r);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
