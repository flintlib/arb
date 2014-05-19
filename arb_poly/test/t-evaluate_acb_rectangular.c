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

    printf("evaluate_acb_rectangular....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        arb_poly_t f;
        acb_t x, y1, y2;

        acb_init(x);
        acb_init(y1);
        acb_init(y2);
        arb_poly_init(f);

        acb_randtest(x, state, 2 + n_randint(state, 1000), 5);
        arb_poly_randtest(f, state, 2 + n_randint(state, 100), 2 + n_randint(state, 1000), 5);

        arb_poly_evaluate_acb_rectangular(y1, f, x, 2 + n_randint(state, 1000));
        arb_poly_evaluate_acb_horner(y2, f, x, 2 + n_randint(state, 1000));

        if (!acb_overlaps(y1, y2))
        {
            printf("FAIL\n\n");
            printf("f = "); arb_poly_printd(f, 15); printf("\n\n");
            printf("x = "); acb_printd(x, 15); printf("\n\n");
            printf("y1 = "); acb_printd(y1, 15); printf("\n\n");
            printf("y2 = "); acb_printd(y2, 15); printf("\n\n");
            abort();
        }

        arb_poly_clear(f);
        acb_clear(x);
        acb_clear(y1);
        acb_clear(y2);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

