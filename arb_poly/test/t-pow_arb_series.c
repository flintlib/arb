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

#include "arb_poly.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("pow_arb_series....");
    fflush(stdout);

    flint_randinit(state);

    /* compare with exp/log */
    for (iter = 0; iter < 10000; iter++)
    {
        long prec, trunc;
        arb_poly_t f, g, h1, h2;
        arb_t c;

        prec = 2 + n_randint(state, 200);
        trunc = n_randint(state, 20);

        arb_poly_init(f);
        arb_poly_init(g);
        arb_poly_init(h1);
        arb_poly_init(h2);
        arb_init(c);

        /* generate binomials */
        if (n_randint(state, 20) == 0)
        {
            arb_randtest(c, state, prec, 10);
            arb_poly_set_coeff_arb(f, 0, c);
            arb_randtest(c, state, prec, 10);
            arb_poly_set_coeff_arb(f, 1 + n_randint(state, 20), c);
        }
        else
        {
            arb_poly_randtest(f, state, 1 + n_randint(state, 20), prec, 10);
        }

        arb_poly_randtest(h1, state, 1 + n_randint(state, 20), prec, 10);

        arb_randtest(c, state, prec, 10);
        arb_poly_set_arb(g, c);

        /* f^c */
        arb_poly_pow_arb_series(h1, f, c, trunc, prec);

        /* f^c = exp(c*log(f)) */
        arb_poly_log_series(h2, f, trunc, prec);
        arb_poly_mullow(h2, h2, g, trunc, prec);
        arb_poly_exp_series(h2, h2, trunc, prec);

        if (!arb_poly_overlaps(h1, h2))
        {
            printf("FAIL\n\n");

            printf("prec = %ld\n", prec);
            printf("trunc = %ld\n", trunc);

            printf("f = "); arb_poly_printd(f, 15); printf("\n\n");
            printf("c = "); arb_printd(c, 15); printf("\n\n");
            printf("h1 = "); arb_poly_printd(h1, 15); printf("\n\n");
            printf("h2 = "); arb_poly_printd(h2, 15); printf("\n\n");

            abort();
        }

        arb_poly_pow_arb_series(f, f, c, trunc, prec);

        if (!arb_poly_overlaps(f, h1))
        {
            printf("FAIL (aliasing)\n\n");
            abort();
        }

        arb_poly_clear(f);
        arb_poly_clear(g);
        arb_poly_clear(h1);
        arb_poly_clear(h2);
        arb_clear(c);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

