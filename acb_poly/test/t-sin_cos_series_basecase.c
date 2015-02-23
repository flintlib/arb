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

#include "acb_poly.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("sin_cos_series_basecase....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        long m, n, rbits1, rbits2;
        fmpq_poly_t B;
        acb_poly_t a, b, c, d, e;
        int times_pi;

        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);

        m = 1 + n_randint(state, 30);
        n = 1 + n_randint(state, 30);
        times_pi = n_randint(state, 2);

        fmpq_poly_init(B);
        acb_poly_init(a);
        acb_poly_init(b);
        acb_poly_init(c);
        acb_poly_init(d);
        acb_poly_init(e);

        acb_poly_randtest(a, state, m, rbits1, 10);

        acb_poly_sin_cos_series_basecase(b, c, a, n, rbits2, times_pi);

        /* Check sin(x)^2 + cos(x)^2 = 1 */
        acb_poly_mullow(d, b, b, n, rbits2);
        acb_poly_mullow(e, c, c, n, rbits2);
        acb_poly_add(d, d, e, rbits2);

        fmpq_poly_one(B);
        if (!acb_poly_contains_fmpq_poly(d, B))
        {
            printf("FAIL\n\n");
            printf("bits2 = %ld\n", rbits2);

            printf("a = "); acb_poly_printd(a, 15); printf("\n\n");
            printf("b = "); acb_poly_printd(b, 15); printf("\n\n");
            printf("c = "); acb_poly_printd(c, 15); printf("\n\n");
            printf("d = "); acb_poly_printd(d, 15); printf("\n\n");

            abort();
        }

        acb_poly_set(d, a);
        acb_poly_sin_cos_series_basecase(d, c, d, n, rbits2, times_pi);
        if (!acb_poly_equal(b, d))
        {
            printf("FAIL (aliasing 1)\n\n");
            abort();
        }

        acb_poly_set(d, a);
        acb_poly_sin_cos_series_basecase(b, d, d, n, rbits2, times_pi);
        if (!acb_poly_equal(c, d))
        {
            printf("FAIL (aliasing 2)\n\n");
            abort();
        }

        fmpq_poly_clear(B);
        acb_poly_clear(a);
        acb_poly_clear(b);
        acb_poly_clear(c);
        acb_poly_clear(d);
        acb_poly_clear(e);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

