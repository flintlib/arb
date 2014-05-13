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

    printf("pow_series....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 50000; iter++)
    {
        long rbits1, rbits2, len;
        arb_poly_t a, b, c, d;

        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);
        len = n_randint(state, 25);

        arb_poly_init(a);
        arb_poly_init(b);
        arb_poly_init(c);
        arb_poly_init(d);

        if (n_randint(state, 4) == 0)
            arb_poly_randtest(a, state, 1, rbits1, 25);
        else
            arb_poly_randtest(a, state, 1 + n_randint(state, 20), rbits1, 5);

        if (n_randint(state, 4) == 0)
            arb_poly_randtest(b, state, 1, rbits1, 25);
        else
            arb_poly_randtest(b, state, 1 + n_randint(state, 20), rbits1, 5);

        arb_poly_randtest(c, state, 1 + n_randint(state, 20), rbits1, 5);

        arb_poly_pow_series(c, a, b, len, rbits2);

        /* a^b = exp(b log a) */
        arb_poly_log_series(d, a, len, rbits2);
        arb_poly_mullow(d, d, b, len, rbits2);
        arb_poly_exp_series(d, d, len, rbits2);

        if (!arb_poly_overlaps(c, d))
        {
            printf("FAIL (iter %ld)\n\n", iter);
            printf("bits2 = %ld\n", rbits2);
            printf("len = %ld\n", len);

            printf("a = "); arb_poly_printd(a, 15); printf("\n\n");
            printf("b = "); arb_poly_printd(b, 15); printf("\n\n");
            printf("c = "); arb_poly_printd(c, 15); printf("\n\n");
            printf("d = "); arb_poly_printd(d, 15); printf("\n\n");

            abort();
        }

        /* check aliasing */
        if (iter < 5000)
        {
            arb_poly_set(d, a);
            arb_poly_pow_series(d, d, b, len, rbits2);

            if (!arb_poly_overlaps(c, d))
            {
                printf("FAIL (aliasing 1)\n");
                abort();
            }

            arb_poly_set(d, b);
            arb_poly_pow_series(d, a, d, len, rbits2);

            if (!arb_poly_overlaps(c, d))
            {
                printf("FAIL (aliasing 2)\n");
                abort();
            }
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
