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

#include "acb_poly.h"


int main()
{
    long iter;
    flint_rand_t state;

    printf("find_roots....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        acb_poly_t A;
        acb_poly_t B;
        acb_poly_t C;
        acb_t t;
        acb_ptr roots;
        long i, deg, isolated;
        long prec = 10 + n_randint(state, 400);

        acb_init(t);
        acb_poly_init(A);
        acb_poly_init(B);
        acb_poly_init(C);

        do {
            acb_poly_randtest(A, state, 2 + n_randint(state, 15), prec, 5);
        } while (A->length == 0);
        deg = A->length - 1;

        roots = _acb_vec_init(deg);

        isolated = acb_poly_find_roots(roots, A, NULL, 0, prec);

        if (isolated == deg)
        {
            acb_poly_fit_length(B, 1);
            acb_set(B->coeffs, A->coeffs + deg);
            _acb_poly_set_length(B, 1);

            for (i = 0; i < deg; i++)
            {
                acb_poly_fit_length(C, 2);
                acb_one(C->coeffs + 1);
                acb_neg(C->coeffs + 0, roots + i);
                _acb_poly_set_length(C, 2);
                acb_poly_mul(B, B, C, prec);
            }

            if (!acb_poly_contains(B, A))
            {
                printf("FAIL: product does not equal polynomial\n");
                acb_poly_printd(A, 15); printf("\n\n");
                acb_poly_printd(B, 15); printf("\n\n");
                abort();
            }
        }

        for (i = 0; i < isolated; i++)
        {
            acb_poly_evaluate(t, A, roots + i, prec);
            if (!acb_contains_zero(t))
            {
                printf("FAIL: poly(root) does not contain zero\n");
                acb_poly_printd(A, 15); printf("\n\n");
                acb_printd(roots + i, 15); printf("\n\n");
                acb_printd(t, 15); printf("\n\n");
                abort();
            }
        }

        _acb_vec_clear(roots, deg);

        acb_clear(t);
        acb_poly_clear(A);
        acb_poly_clear(B);
        acb_poly_clear(C);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
