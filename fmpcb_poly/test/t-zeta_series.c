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

    printf("zeta_series....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 400; iter++)
    {
        long m, n1, n2, bits1, bits2, bits3;
        int deflate;
        fmpcb_poly_t S, A, B, C, D, E, F;
        fmpcb_t a, a1;

        bits1 = 2 + n_randint(state, 200);
        bits2 = 2 + n_randint(state, 200);
        bits3 = 2 + n_randint(state, 200);

        m = 1 + n_randint(state, 25);
        n1 = 1 + n_randint(state, 25);
        n2 = 1 + n_randint(state, 25);

        fmpcb_poly_init(S);
        fmpcb_poly_init(A);
        fmpcb_poly_init(B);
        fmpcb_poly_init(C);
        fmpcb_poly_init(D);
        fmpcb_poly_init(E);
        fmpcb_poly_init(F);
        fmpcb_init(a);
        fmpcb_init(a1);

        deflate = n_randint(state, 2);

        fmpcb_poly_randtest(S, state, m, bits1, 3);
        fmprb_randtest_precise(fmpcb_realref(a), state, bits1, 3);
        fmprb_randtest_precise(fmpcb_imagref(a), state, bits1, 3);
        fmpcb_poly_set_coeff_fmpcb(S, 0, a);

        if (n_randint(state, 2))
            fmpcb_randtest(a, state, bits1, 3);
        else
            fmpcb_one(a);

        fmpcb_poly_zeta_series(A, S, a, deflate, n1, bits2);
        fmpcb_poly_zeta_series(B, S, a, deflate, n2, bits3);

        fmpcb_poly_set(C, A);
        fmpcb_poly_truncate(C, FLINT_MIN(n1, n2));
        fmpcb_poly_truncate(B, FLINT_MIN(n1, n2));

        if (!fmpcb_poly_overlaps(B, C))
        {
            printf("FAIL\n\n");
            printf("S = "); fmpcb_poly_printd(S, 15); printf("\n\n");
            printf("a = "); fmpcb_printd(a, 15); printf("\n\n");
            printf("A = "); fmpcb_poly_printd(A, 15); printf("\n\n");
            printf("B = "); fmpcb_poly_printd(B, 15); printf("\n\n");
            abort();
        }

        /* check zeta(s,a) = zeta(s,a+1) + a^(-s) */
        fmpcb_poly_set_fmpcb(D, a);
        fmpcb_poly_log_series(D, D, n1, bits2);
        fmpcb_poly_mullow(D, D, S, n1, bits2);
        fmpcb_poly_neg(D, D);
        fmpcb_poly_exp_series(D, D, n1, bits2);

        fmpcb_add_ui(a1, a, 1, bits2);
        fmpcb_poly_zeta_series(E, S, a1, deflate, n1, bits2);
        fmpcb_poly_add(E, E, D, bits2);

        if (!fmpcb_poly_overlaps(A, E))
        {
            printf("FAIL (functional equation)\n\n");
            printf("S = "); fmpcb_poly_printd(S, 15); printf("\n\n");
            printf("a = "); fmpcb_printd(a, 15); printf("\n\n");
            printf("A = "); fmpcb_poly_printd(A, 15); printf("\n\n");
            printf("E = "); fmpcb_poly_printd(A, 15); printf("\n\n");
            abort();
        }

        fmpcb_poly_zeta_series(S, S, a, deflate, n1, bits2);
        if (!fmpcb_poly_overlaps(A, S))
        {
            printf("FAIL (aliasing)\n\n");
            abort();
        }

        fmpcb_poly_clear(S);
        fmpcb_poly_clear(A);
        fmpcb_poly_clear(B);
        fmpcb_poly_clear(C);
        fmpcb_poly_clear(D);
        fmpcb_poly_clear(E);
        fmpcb_poly_clear(F);
        fmpcb_clear(a);
        fmpcb_clear(a1);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

