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

#include "acb_poly.h"


int main()
{
    long iter;
    flint_rand_t state;

    printf("mullow_transpose....");
    fflush(stdout);

    flint_randinit(state);

    /* compare with fmpq_poly */
    for (iter = 0; iter < 10000; iter++)
    {
        long qbits1, qbits2, rbits1, rbits2, rbits3, trunc;
        fmpq_poly_t A, B, C;
        acb_poly_t a, b, c, d;

        qbits1 = 2 + n_randint(state, 200);
        qbits2 = 2 + n_randint(state, 200);
        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);
        rbits3 = 2 + n_randint(state, 200);
        trunc = n_randint(state, 10);

        fmpq_poly_init(A);
        fmpq_poly_init(B);
        fmpq_poly_init(C);

        acb_poly_init(a);
        acb_poly_init(b);
        acb_poly_init(c);
        acb_poly_init(d);

        fmpq_poly_randtest(A, state, 1 + n_randint(state, 10), qbits1);
        fmpq_poly_randtest(B, state, 1 + n_randint(state, 10), qbits2);

        fmpq_poly_mullow(C, A, B, trunc);

        acb_poly_set_fmpq_poly(a, A, rbits1);
        acb_poly_set_fmpq_poly(b, B, rbits2);

        acb_poly_mullow_transpose(c, a, b, trunc, rbits3);

        if (!acb_poly_contains_fmpq_poly(c, C))
        {
            printf("FAIL\n\n");
            printf("bits3 = %ld\n", rbits3);
            printf("trunc = %ld\n", trunc);

            printf("A = "); fmpq_poly_print(A); printf("\n\n");
            printf("B = "); fmpq_poly_print(B); printf("\n\n");
            printf("C = "); fmpq_poly_print(C); printf("\n\n");

            printf("a = "); acb_poly_printd(a, 15); printf("\n\n");
            printf("b = "); acb_poly_printd(b, 15); printf("\n\n");
            printf("c = "); acb_poly_printd(c, 15); printf("\n\n");

            abort();
        }

        acb_poly_set(d, a);
        acb_poly_mullow_transpose(d, d, b, trunc, rbits3);
        if (!acb_poly_equal(d, c))
        {
            printf("FAIL (aliasing 1)\n\n");
            abort();
        }

        acb_poly_set(d, b);
        acb_poly_mullow_transpose(d, a, d, trunc, rbits3);
        if (!acb_poly_equal(d, c))
        {
            printf("FAIL (aliasing 2)\n\n");
            abort();
        }

        /* test squaring */
        acb_poly_set(b, a);
        acb_poly_mullow_transpose(c, a, b, trunc, rbits3);
        acb_poly_mullow_transpose(d, a, a, trunc, rbits3);
        if (!acb_poly_overlaps(c, d))  /* not guaranteed to be identical */
        {
            printf("FAIL (squaring)\n\n");

            printf("a = "); acb_poly_printd(a, 15); printf("\n\n");
            printf("b = "); acb_poly_printd(b, 15); printf("\n\n");
            printf("c = "); acb_poly_printd(c, 15); printf("\n\n");

            abort();
        }

        acb_poly_mullow_transpose(a, a, a, trunc, rbits3);
        if (!acb_poly_equal(d, a))
        {
            printf("FAIL (aliasing, squaring)\n\n");

            printf("a = "); acb_poly_printd(a, 15); printf("\n\n");
            printf("b = "); acb_poly_printd(b, 15); printf("\n\n");
            printf("d = "); acb_poly_printd(d, 15); printf("\n\n");

            abort();
        }

        fmpq_poly_clear(A);
        fmpq_poly_clear(B);
        fmpq_poly_clear(C);

        acb_poly_clear(a);
        acb_poly_clear(b);
        acb_poly_clear(c);
        acb_poly_clear(d);
    }

    /* compare with classical */
    for (iter = 0; iter < 10000; iter++)
    {
        long bits, trunc;
        acb_poly_t a, b, ab, ab2;

        bits = 2 + n_randint(state, 200);
        trunc = n_randint(state, 10);

        acb_poly_init(a);
        acb_poly_init(b);
        acb_poly_init(ab);
        acb_poly_init(ab2);

        acb_poly_randtest(a, state, 1 + n_randint(state, 10), bits, 5);
        acb_poly_randtest(b, state, 1 + n_randint(state, 10), bits, 5);

        acb_poly_mullow_classical(ab, a, b, trunc, bits);
        acb_poly_mullow_transpose(ab2, a, b, trunc, bits);

        if (!acb_poly_overlaps(ab, ab2))
        {
            printf("FAIL\n\n");
            printf("bits = %ld\n", bits);
            printf("trunc = %ld\n", trunc);

            printf("a = "); acb_poly_printd(a, 15); printf("\n\n");
            printf("b = "); acb_poly_printd(b, 15); printf("\n\n");
            printf("ab = "); acb_poly_printd(ab, 15); printf("\n\n");
            printf("ab2 = "); acb_poly_printd(ab2, 15); printf("\n\n");

            abort();
        }

        acb_poly_clear(a);
        acb_poly_clear(b);
        acb_poly_clear(ab);
        acb_poly_clear(ab2);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

