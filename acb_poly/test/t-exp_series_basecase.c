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

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include "acb_poly.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("exp_series_basecase....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        long m, n, bits1, bits2, bits3;
        acb_poly_t a, b, c, d;

        bits1 = 2 + n_randint(state, 200);
        bits2 = 2 + n_randint(state, 200);
        bits3 = 2 + n_randint(state, 200);

        m = 1 + n_randint(state, 40);
        n = 1 + n_randint(state, 40);

        acb_poly_init(a);
        acb_poly_init(b);
        acb_poly_init(c);
        acb_poly_init(d);

        acb_poly_randtest(a, state, m, bits1, 10);
        acb_poly_randtest(b, state, m, bits1, 10);

        /* check exp(a+b) = exp(a) exp(b) */
        acb_poly_exp_series_basecase(c, a, n, bits2);
        acb_poly_exp_series_basecase(d, b, n, bits2);
        acb_poly_mullow(c, c, d, n, bits2);

        acb_poly_add(d, a, b, bits3);
        acb_poly_exp_series_basecase(d, d, n, bits3);   /* also aliasing test */

        if (!acb_poly_overlaps(c, d))
        {
            printf("FAIL\n\n");

            printf("a = "); acb_poly_printd(a, 15); printf("\n\n");
            printf("b = "); acb_poly_printd(b, 15); printf("\n\n");
            printf("c = "); acb_poly_printd(c, 15); printf("\n\n");
            printf("d = "); acb_poly_printd(d, 15); printf("\n\n");
            abort();
        }

        acb_poly_clear(a);
        acb_poly_clear(b);
        acb_poly_clear(c);
        acb_poly_clear(d);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
