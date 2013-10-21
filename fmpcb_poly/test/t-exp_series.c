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

#include "fmpcb_poly.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("exp_series....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        long m, n, bits1, bits2, bits3;
        fmpcb_poly_t a, b, c, d;

        bits1 = 2 + n_randint(state, 200);
        bits2 = 2 + n_randint(state, 200);
        bits3 = 2 + n_randint(state, 200);

        if (n_randint(state, 100) == 0)
        {
            m = 1 + n_randint(state, 300);
            n = 1 + n_randint(state, 300);
        }
        else
        {
            m = 1 + n_randint(state, 20);
            n = 1 + n_randint(state, 20);
        }

        fmpcb_poly_init(a);
        fmpcb_poly_init(b);
        fmpcb_poly_init(c);
        fmpcb_poly_init(d);

        fmpcb_poly_randtest(a, state, m, bits1, 10);
        fmpcb_poly_randtest(b, state, m, bits1, 10);

        fmpcb_poly_randtest(c, state, 1 + n_randint(state, 300), bits1, 10);
        fmpcb_poly_randtest(d, state, 1 + n_randint(state, 300), bits1, 10);

        /* check exp(a+b) = exp(a) exp(b) */
        fmpcb_poly_exp_series(c, a, n, bits2);
        fmpcb_poly_exp_series(d, b, n, bits2);
        fmpcb_poly_mullow(c, c, d, n, bits2);

        fmpcb_poly_add(d, a, b, bits3);
        fmpcb_poly_exp_series(d, d, n, bits3);   /* also aliasing test */

        if (!fmpcb_poly_overlaps(c, d))
        {
            printf("FAIL\n\n");

            printf("a = "); fmpcb_poly_printd(a, 15); printf("\n\n");
            printf("b = "); fmpcb_poly_printd(b, 15); printf("\n\n");
            printf("c = "); fmpcb_poly_printd(c, 15); printf("\n\n");
            printf("d = "); fmpcb_poly_printd(d, 15); printf("\n\n");
            abort();
        }

        fmpcb_poly_clear(a);
        fmpcb_poly_clear(b);
        fmpcb_poly_clear(c);
        fmpcb_poly_clear(d);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
