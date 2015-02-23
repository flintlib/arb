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

    printf("sin_pi_series....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        long m, n1, n2, bits1, bits2, bits3;
        arb_poly_t S, A, B, C;
        arb_t pi;

        bits1 = 2 + n_randint(state, 200);
        bits2 = 2 + n_randint(state, 200);
        bits3 = 2 + n_randint(state, 200);

        m = 1 + n_randint(state, 30);
        n1 = 1 + n_randint(state, 30);
        n2 = 1 + n_randint(state, 30);

        arb_poly_init(S);
        arb_poly_init(A);
        arb_poly_init(B);
        arb_poly_init(C);
        arb_init(pi);

        arb_poly_randtest(S, state, m, bits1, 3);
        arb_poly_randtest(A, state, m, bits1, 3);
        arb_poly_randtest(B, state, m, bits1, 3);

        arb_poly_sin_pi_series(A, S, n1, bits2);

        arb_const_pi(pi, bits3);
        arb_poly_set_arb(B, pi);
        arb_poly_mul(B, S, B, bits3);
        arb_poly_sin_series(B, B, n2, bits3);

        arb_poly_set(C, A);
        arb_poly_truncate(C, FLINT_MIN(n1, n2));
        arb_poly_truncate(B, FLINT_MIN(n1, n2));

        if (!arb_poly_overlaps(B, C))
        {
            printf("FAIL\n\n");
            printf("S = "); arb_poly_printd(S, 15); printf("\n\n");
            printf("A = "); arb_poly_printd(A, 15); printf("\n\n");
            printf("B = "); arb_poly_printd(B, 15); printf("\n\n");
            abort();
        }

        arb_poly_sin_pi_series(S, S, n1, bits2);

        if (!arb_poly_overlaps(A, S))
        {
            printf("FAIL (aliasing)\n\n");
            abort();
        }

        arb_poly_clear(S);
        arb_poly_clear(A);
        arb_poly_clear(B);
        arb_poly_clear(C);
        arb_clear(pi);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
