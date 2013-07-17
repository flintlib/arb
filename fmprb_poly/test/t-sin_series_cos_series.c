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

#include "fmprb_poly.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("sin_series/cos_series....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 2000; iter++)
    {
        long m, n, qbits, rbits1, rbits2;
        fmpq_poly_t A, B;
        fmprb_poly_t a, b, c, d, e;

        qbits = 2 + n_randint(state, 200);
        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);

        m = 1 + n_randint(state, 30);
        n = 1 + n_randint(state, 30);

        fmpq_poly_init(A);
        fmpq_poly_init(B);
        fmprb_poly_init(a);
        fmprb_poly_init(b);
        fmprb_poly_init(c);
        fmprb_poly_init(d);
        fmprb_poly_init(e);

        fmpq_poly_randtest(A, state, m, qbits);
        fmprb_poly_set_fmpq_poly(a, A, rbits1);

        fmprb_poly_sin_series(b, a, n, rbits2);
        fmprb_poly_cos_series(c, a, n, rbits2);

        /* Check sin(x)^2 + cos(x)^2 = 1 */
        fmprb_poly_mullow(d, b, b, n, rbits2);
        fmprb_poly_mullow(e, c, c, n, rbits2);
        fmprb_poly_add(d, d, e, rbits2);

        fmpq_poly_one(B);
        if (!fmprb_poly_contains_fmpq_poly(d, B))
        {
            printf("FAIL\n\n");
            printf("bits2 = %ld\n", rbits2);

            printf("A = "); fmpq_poly_print(A); printf("\n\n");
            printf("a = "); fmprb_poly_printd(a, 15); printf("\n\n");
            printf("b = "); fmprb_poly_printd(b, 15); printf("\n\n");
            printf("c = "); fmprb_poly_printd(c, 15); printf("\n\n");
            printf("d = "); fmprb_poly_printd(d, 15); printf("\n\n");

            abort();
        }

        fmprb_poly_set(d, a);
        fmprb_poly_sin_series(d, d, n, rbits2);
        if (!fmprb_poly_equal(b, d))
        {
            printf("FAIL (aliasing 1)\n\n");
            abort();
        }

        fmprb_poly_set(d, a);
        fmprb_poly_cos_series(d, d, n, rbits2);
        if (!fmprb_poly_equal(c, d))
        {
            printf("FAIL (aliasing 2)\n\n");
            abort();
        }

        fmpq_poly_clear(A);
        fmpq_poly_clear(B);
        fmprb_poly_clear(a);
        fmprb_poly_clear(b);
        fmprb_poly_clear(c);
        fmprb_poly_clear(d);
        fmprb_poly_clear(e);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

