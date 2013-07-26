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

#include "fmprb_poly.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("zeta_series....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 500; iter++)
    {
        long m, n1, n2, bits1, bits2, bits3;
        int deflate;
        fmprb_poly_t S, A, B, C, D, E, F;
        fmprb_t a, a1;

        bits1 = 2 + n_randint(state, 300);
        bits2 = 2 + n_randint(state, 300);
        bits3 = 2 + n_randint(state, 300);

        m = 1 + n_randint(state, 30);
        n1 = 1 + n_randint(state, 30);
        n2 = 1 + n_randint(state, 30);

        fmprb_poly_init(S);
        fmprb_poly_init(A);
        fmprb_poly_init(B);
        fmprb_poly_init(C);
        fmprb_poly_init(D);
        fmprb_poly_init(E);
        fmprb_poly_init(F);
        fmprb_init(a);
        fmprb_init(a1);

        deflate = n_randint(state, 2);

        fmprb_poly_randtest(S, state, m, bits1, 3);
        fmprb_randtest_precise(a, state, bits1, 3);
        fmprb_poly_set_coeff_fmprb(S, 0, a);

        if (n_randint(state, 2))
            fmprb_randtest(a, state, bits1, 3);
        else
            fmprb_one(a);

        fmprb_poly_zeta_series(A, S, a, deflate, n1, bits2);
        fmprb_poly_zeta_series(B, S, a, deflate, n2, bits3);

        fmprb_poly_set(C, A);
        fmprb_poly_truncate(C, FLINT_MIN(n1, n2));
        fmprb_poly_truncate(B, FLINT_MIN(n1, n2));

        if (!fmprb_poly_overlaps(B, C))
        {
            printf("FAIL\n\n");
            printf("S = "); fmprb_poly_printd(S, 15); printf("\n\n");
            printf("a = "); fmprb_printd(a, 15); printf("\n\n");
            printf("A = "); fmprb_poly_printd(A, 15); printf("\n\n");
            printf("B = "); fmprb_poly_printd(B, 15); printf("\n\n");
            abort();
        }

        /* check zeta(s,a) = zeta(s,a+1) + a^(-s) */
        fmprb_poly_set_fmprb(D, a);
        fmprb_poly_log_series(D, D, n1, bits2);
        fmprb_poly_mullow(D, D, S, n1, bits2);
        fmprb_poly_neg(D, D);
        fmprb_poly_exp_series(D, D, n1, bits2);

        fmprb_add_ui(a1, a, 1, bits2);
        fmprb_poly_zeta_series(E, S, a1, deflate, n1, bits2);
        fmprb_poly_add(E, E, D, bits2);

        if (!fmprb_poly_overlaps(A, E))
        {
            printf("FAIL (functional equation)\n\n");
            printf("S = "); fmprb_poly_printd(S, 15); printf("\n\n");
            printf("a = "); fmprb_printd(a, 15); printf("\n\n");
            printf("A = "); fmprb_poly_printd(A, 15); printf("\n\n");
            printf("E = "); fmprb_poly_printd(A, 15); printf("\n\n");
            abort();
        }

        fmprb_poly_zeta_series(S, S, a, deflate, n1, bits2);
        if (!fmprb_poly_overlaps(A, S))
        {
            printf("FAIL (aliasing)\n\n");
            abort();
        }

        fmprb_poly_clear(S);
        fmprb_poly_clear(A);
        fmprb_poly_clear(B);
        fmprb_poly_clear(C);
        fmprb_poly_clear(D);
        fmprb_poly_clear(E);
        fmprb_poly_clear(F);
        fmprb_clear(a);
        fmprb_clear(a1);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

