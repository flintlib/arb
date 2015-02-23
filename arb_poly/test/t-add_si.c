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

#include "arb_poly.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("add_si....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        arb_poly_t a, b, c, d;
        long v;

        arb_poly_init(a);
        arb_poly_init(b);
        arb_poly_init(c);
        arb_poly_init(d);

        arb_poly_randtest(a, state, 1 + n_randint(state, 10), 1 + n_randint(state, 200), 10);
        arb_poly_randtest(b, state, 1 + n_randint(state, 10), 1 + n_randint(state, 200), 10);
        arb_poly_randtest(c, state, 1 + n_randint(state, 10), 1 + n_randint(state, 200), 10);
        arb_poly_randtest(d, state, 1 + n_randint(state, 10), 1 + n_randint(state, 200), 10);

        v = n_randtest(state);

        arb_poly_set_si(b, v);

        arb_poly_add(c, a, b, 2 + n_randint(state, 200));
        arb_poly_add_si(d, a, v, 2 + n_randint(state, 200));

        if (!arb_poly_overlaps(c, d))
        {
            printf("FAIL\n\n");
            abort();
        }

        arb_poly_add_si(a, a, v, 2 + n_randint(state, 200));

        if (!arb_poly_overlaps(a, d))
        {
            printf("FAIL (aliasing)\n\n");
            abort();
        }

        arb_poly_clear(a);
        arb_poly_clear(b);
        arb_poly_clear(c);
        arb_poly_clear(d);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

