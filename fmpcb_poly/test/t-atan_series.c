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

#include "fmpcb_poly.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("atan_series....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        long m, n, qbits, rbits1, rbits2;
        fmpq_poly_t A;
        fmpcb_poly_t a, b, c, d;

        qbits = 2 + n_randint(state, 200);
        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);

        m = 1 + n_randint(state, 30);
        n = 1 + n_randint(state, 30);

        fmpq_poly_init(A);
        fmpcb_poly_init(a);
        fmpcb_poly_init(b);
        fmpcb_poly_init(c);
        fmpcb_poly_init(d);

        fmpq_poly_randtest(A, state, m, qbits);
        fmpq_poly_set_coeff_si(A, 0, 0);
        fmpcb_poly_set_fmpq_poly(a, A, rbits1);

        fmpcb_poly_atan_series(b, a, n, rbits2);

        /* Check 2 atan(x) = atan(2x/(1-x^2)) + C */
        fmpcb_poly_mullow(c, a, a, n, rbits2);
        fmpcb_poly_one(d);
        fmpcb_poly_sub(c, d, c, rbits2);
        fmpcb_poly_add(d, a, a, rbits2);

        if (fmpcb_poly_length(c) != 0)
        {
            fmpcb_poly_div_series(c, d, c, n, rbits2);
            fmpcb_poly_atan_series(c, c, n, rbits2);
            fmpcb_poly_add(d, b, b, rbits2);

            /* TODO: also check the first coefficient */
            fmpcb_poly_set_coeff_si(c, 0, 0);
            fmpcb_poly_set_coeff_si(d, 0, 0);

            if (!fmpcb_poly_overlaps(c, d))
            {
                printf("FAIL\n\n");
                printf("bits2 = %ld\n", rbits2);

                printf("A = "); fmpq_poly_print(A); printf("\n\n");
                printf("a = "); fmpcb_poly_printd(a, 15); printf("\n\n");
                printf("b = "); fmpcb_poly_printd(b, 15); printf("\n\n");
                printf("c = "); fmpcb_poly_printd(c, 15); printf("\n\n");
                printf("d = "); fmpcb_poly_printd(d, 15); printf("\n\n");

                abort();
            }
        }

        fmpcb_poly_atan_series(a, a, n, rbits2);
        if (!fmpcb_poly_equal(a, b))
        {
            printf("FAIL (aliasing)\n\n");
            abort();
        }

        fmpq_poly_clear(A);
        fmpcb_poly_clear(a);
        fmpcb_poly_clear(b);
        fmpcb_poly_clear(c);
        fmpcb_poly_clear(d);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

