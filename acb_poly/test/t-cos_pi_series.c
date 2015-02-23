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

    printf("cos_pi_series....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        long m, n1, n2, bits1, bits2, bits3;
        acb_poly_t S, A, B, C;
        acb_t pi;

        bits1 = 2 + n_randint(state, 200);
        bits2 = 2 + n_randint(state, 200);
        bits3 = 2 + n_randint(state, 200);

        m = 1 + n_randint(state, 30);
        n1 = 1 + n_randint(state, 30);
        n2 = 1 + n_randint(state, 30);

        acb_poly_init(S);
        acb_poly_init(A);
        acb_poly_init(B);
        acb_poly_init(C);
        acb_init(pi);

        acb_poly_randtest(S, state, m, bits1, 3);
        acb_poly_randtest(A, state, m, bits1, 3);
        acb_poly_randtest(B, state, m, bits1, 3);

        acb_poly_cos_pi_series(A, S, n1, bits2);

        acb_const_pi(pi, bits3);
        acb_poly_set_acb(B, pi);
        acb_poly_mul(B, S, B, bits3);
        acb_poly_cos_series(B, B, n2, bits3);

        acb_poly_set(C, A);
        acb_poly_truncate(C, FLINT_MIN(n1, n2));
        acb_poly_truncate(B, FLINT_MIN(n1, n2));

        if (!acb_poly_overlaps(B, C))
        {
            printf("FAIL\n\n");
            printf("S = "); acb_poly_printd(S, 15); printf("\n\n");
            printf("A = "); acb_poly_printd(A, 15); printf("\n\n");
            printf("B = "); acb_poly_printd(B, 15); printf("\n\n");
            abort();
        }

        acb_poly_cos_pi_series(S, S, n1, bits2);

        if (!acb_poly_overlaps(A, S))
        {
            printf("FAIL (aliasing)\n\n");
            abort();
        }

        acb_poly_clear(S);
        acb_poly_clear(A);
        acb_poly_clear(B);
        acb_poly_clear(C);
        acb_clear(pi);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
