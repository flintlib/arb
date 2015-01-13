/*=============================================================================

    This file is part of acb.

    acb is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    acb is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with acb; if not, write to the Free Software
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

    printf("pow_acb_series....");
    fflush(stdout);

    flint_randinit(state);

    /* compare with exp/log */
    for (iter = 0; iter < 10000; iter++)
    {
        long prec, trunc;
        acb_poly_t f, g, h1, h2;
        acb_t c;

        prec = 2 + n_randint(state, 200);
        trunc = n_randint(state, 20);

        acb_poly_init(f);
        acb_poly_init(g);
        acb_poly_init(h1);
        acb_poly_init(h2);
        acb_init(c);

        /* generate binomials */
        if (n_randint(state, 20) == 0)
        {
            acb_randtest(c, state, prec, 10);
            acb_poly_set_coeff_acb(f, 0, c);
            acb_randtest(c, state, prec, 10);
            acb_poly_set_coeff_acb(f, 1 + n_randint(state, 20), c);
        }
        else
        {
            acb_poly_randtest(f, state, 1 + n_randint(state, 20), prec, 10);
        }

        acb_poly_randtest(h1, state, 1 + n_randint(state, 20), prec, 10);

        acb_randtest(c, state, prec, 10);
        acb_poly_set_acb(g, c);

        /* f^c */
        acb_poly_pow_acb_series(h1, f, c, trunc, prec);

        /* f^c = exp(c*log(f)) */
        acb_poly_log_series(h2, f, trunc, prec);
        acb_poly_mullow(h2, h2, g, trunc, prec);
        acb_poly_exp_series(h2, h2, trunc, prec);

        if (!acb_poly_overlaps(h1, h2))
        {
            printf("FAIL\n\n");

            printf("prec = %ld\n", prec);
            printf("trunc = %ld\n", trunc);

            printf("f = "); acb_poly_printd(f, 15); printf("\n\n");
            printf("c = "); acb_printd(c, 15); printf("\n\n");
            printf("h1 = "); acb_poly_printd(h1, 15); printf("\n\n");
            printf("h2 = "); acb_poly_printd(h2, 15); printf("\n\n");

            abort();
        }

        acb_poly_pow_acb_series(f, f, c, trunc, prec);

        if (!acb_poly_overlaps(f, h1))
        {
            printf("FAIL (aliasing)\n\n");
            abort();
        }

        acb_poly_clear(f);
        acb_poly_clear(g);
        acb_poly_clear(h1);
        acb_poly_clear(h2);
        acb_clear(c);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

