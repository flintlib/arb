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

#include "fmpcb_poly.h"


int main()
{
    long iter;
    flint_rand_t state;

    printf("rising_ui_series....");
    fflush(stdout);

    flint_randinit(state);

    /* check rf(f, a) * rf(f + a, b) = rf(f, a + b) */
    for (iter = 0; iter < 1000; iter++)
    {
        long bits, trunc;
        ulong a, b;
        fmpcb_poly_t f, g, h1, h2, h1h2, h3;

        bits = 2 + n_randint(state, 200);
        trunc = 1 + n_randint(state, 20);
        a = n_randint(state, 10);
        b = n_randint(state, 10);

        fmpcb_poly_init(f);
        fmpcb_poly_init(g);
        fmpcb_poly_init(h1);
        fmpcb_poly_init(h2);
        fmpcb_poly_init(h1h2);
        fmpcb_poly_init(h3);

        fmpcb_poly_randtest(f, state, 1 + n_randint(state, 20), bits, 4);
        fmpcb_poly_set(g, f);

        /* g = f + 1 */
        if (g->length == 0)
        {
            fmpcb_poly_fit_length(g, 1);
            fmpcb_set_ui(g->coeffs, a);
            _fmpcb_poly_set_length(g, 1);
            _fmpcb_poly_normalise(g);
        }
        else
        {
            fmpcb_add_ui(g->coeffs, g->coeffs, a, bits);
            _fmpcb_poly_normalise(g);
        }

        fmpcb_poly_rising_ui_series(h1, f, a, trunc, bits);
        fmpcb_poly_rising_ui_series(h2, g, b, trunc, bits);
        fmpcb_poly_rising_ui_series(h3, f, a + b, trunc, bits);

        fmpcb_poly_mullow(h1h2, h1, h2, trunc, bits);

        if (!fmpcb_poly_overlaps(h1h2, h3))
        {
            printf("FAIL\n\n");
            printf("bits = %ld\n", bits);
            printf("trunc = %ld\n", trunc);
            printf("a = %lu\n", a);
            printf("b = %lu\n", a);

            printf("f = "); fmpcb_poly_printd(f, 15); printf("\n\n");
            printf("g = "); fmpcb_poly_printd(g, 15); printf("\n\n");
            printf("h1 = "); fmpcb_poly_printd(h1, 15); printf("\n\n");
            printf("h2 = "); fmpcb_poly_printd(h2, 15); printf("\n\n");
            printf("h1h2 = "); fmpcb_poly_printd(h1h2, 15); printf("\n\n");
            printf("h3 = "); fmpcb_poly_printd(h3, 15); printf("\n\n");

            abort();
        }

        fmpcb_poly_rising_ui_series(f, f, a, trunc, bits);

        if (!fmpcb_poly_equal(f, h1))
        {
            printf("FAIL (aliasing)\n\n");

            printf("bits = %ld\n", bits);
            printf("trunc = %ld\n", trunc);
            printf("a = %lu\n", a);
            printf("b = %lu\n", a);

            printf("f = "); fmpcb_poly_printd(f, 15); printf("\n\n");
            printf("h1 = "); fmpcb_poly_printd(h1, 15); printf("\n\n");

            abort();
        }

        fmpcb_poly_clear(f);
        fmpcb_poly_clear(g);
        fmpcb_poly_clear(h1);
        fmpcb_poly_clear(h2);
        fmpcb_poly_clear(h1h2);
        fmpcb_poly_clear(h3);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
