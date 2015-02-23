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

#include "acb_mat.h"

int
main(void)
{
    long iter;
    flint_rand_t state;

    printf("charpoly....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        acb_mat_t A, B, C, D;
        acb_poly_t f, g;
        long m, n;

        m = n_randint(state, 8);
        n = m;

        acb_mat_init(A, m, n);
        acb_mat_init(B, m, n);
        acb_mat_init(C, m, m);
        acb_mat_init(D, n, n);
        acb_poly_init(f);
        acb_poly_init(g);

        acb_mat_randtest(A, state, 1 + n_randint(state, 1000), 10);
        acb_mat_randtest(B, state, 1 + n_randint(state, 1000), 10);

        acb_mat_mul(C, A, B, 2 + n_randint(state, 1000));
        acb_mat_mul(D, B, A, 2 + n_randint(state, 1000));

        acb_mat_charpoly(f, C, 2 + n_randint(state, 1000));
        acb_mat_charpoly(g, D, 2 + n_randint(state, 1000));

        if (!acb_poly_overlaps(f, g))
        {
            flint_printf("FAIL: charpoly(AB) != charpoly(BA).\n");
            flint_printf("Matrix A:\n"), acb_mat_printd(A, 15), flint_printf("\n");
            flint_printf("Matrix B:\n"), acb_mat_printd(B, 15), flint_printf("\n");
            flint_printf("cp(AB) = "), acb_poly_printd(f, 15), flint_printf("\n");
            flint_printf("cp(BA) = "), acb_poly_printd(g, 15), flint_printf("\n");
            abort();
        }

        acb_mat_clear(A);
        acb_mat_clear(B);
        acb_mat_clear(C);
        acb_mat_clear(D);
        acb_poly_clear(f);
        acb_poly_clear(g);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}

