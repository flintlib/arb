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

    printf("tan_series....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 2000; iter++)
    {
        long m, n, qbits, rbits1, rbits2;
        fmpq_poly_t A;
        fmprb_poly_t a, b, c, d, e;

        qbits = 2 + n_randint(state, 200);
        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);

        m = 1 + n_randint(state, 50);
        n = 1 + n_randint(state, 50);

        fmpq_poly_init(A);
        fmprb_poly_init(a);
        fmprb_poly_init(b);
        fmprb_poly_init(c);
        fmprb_poly_init(d);
        fmprb_poly_init(e);

        fmpq_poly_randtest(A, state, m, qbits);
        fmprb_poly_set_fmpq_poly(a, A, rbits1);
        fmprb_poly_randtest(b, state, 1 + n_randint(state, 50), rbits1, 5);

        fmprb_poly_tan_series(b, a, n, rbits2);

        /* check tan(x) = 2*tan(x/2)/(1-tan(x/2)^2) */
        fmprb_poly_scalar_mul_2exp_si(c, a, -1);
        fmprb_poly_tan_series(c, c, n, rbits2);
        fmprb_poly_mullow(d, c, c, n, rbits2);
        fmprb_poly_one(e);
        fmprb_poly_sub(e, e, d, rbits2);
        fmprb_poly_div_series(c, c, e, n, rbits2);
        fmprb_poly_scalar_mul_2exp_si(c, c, 1);

        if (!fmprb_poly_overlaps(b, c))
        {
            printf("FAIL\n\n");
            printf("bits2 = %ld\n", rbits2);

            printf("A = "); fmpq_poly_print(A); printf("\n\n");
            printf("a = "); fmprb_poly_printd(a, 15); printf("\n\n");
            printf("b = "); fmprb_poly_printd(b, 15); printf("\n\n");
            printf("c = "); fmprb_poly_printd(c, 15); printf("\n\n");

            abort();
        }

        fmprb_poly_tan_series(a, a, n, rbits2);
        if (!fmprb_poly_equal(a, b))
        {
            printf("FAIL (aliasing)\n\n");
            abort();
        }

        fmpq_poly_clear(A);
        fmprb_poly_clear(a);
        fmprb_poly_clear(b);
        fmprb_poly_clear(c);
        fmprb_poly_clear(d);
        fmprb_poly_clear(e);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

