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

    printf("evaluate2....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        fmprb_poly_t f, g;
        fmprb_t x, y1, z1, y2, z2;

        fmprb_init(x);
        fmprb_init(y1);
        fmprb_init(z1);
        fmprb_init(y2);
        fmprb_init(z2);
        fmprb_poly_init(f);
        fmprb_poly_init(g);

        fmprb_randtest(x, state, 2 + n_randint(state, 1000), 5);
        fmprb_poly_randtest(f, state, 2 + n_randint(state, 100), 2 + n_randint(state, 1000), 5);
        fmprb_poly_derivative(g, f, 2 + n_randint(state, 1000));

        fmprb_poly_evaluate2(y1, z1, f, x, 2 + n_randint(state, 1000));

        fmprb_poly_evaluate_horner(y2, f, x, 2 + n_randint(state, 1000));
        fmprb_poly_evaluate_horner(z2, g, x, 2 + n_randint(state, 1000));

        if (!fmprb_overlaps(y1, y2) || !fmprb_overlaps(z1, z2))
        {
            printf("FAIL\n\n");
            printf("f = "); fmprb_poly_printd(f, 15); printf("\n\n");
            printf("g = "); fmprb_poly_printd(g, 15); printf("\n\n");
            printf("x = "); fmprb_printd(x, 15); printf("\n\n");
            printf("y1 = "); fmprb_printd(y1, 15); printf("\n\n");
            printf("z1 = "); fmprb_printd(z1, 15); printf("\n\n");
            printf("y2 = "); fmprb_printd(y2, 15); printf("\n\n");
            printf("z2 = "); fmprb_printd(z2, 15); printf("\n\n");
            abort();
        }

        fmprb_poly_clear(f);
        fmprb_poly_clear(g);
        fmprb_clear(x);
        fmprb_clear(y1);
        fmprb_clear(z1);
        fmprb_clear(y2);
        fmprb_clear(z2);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

