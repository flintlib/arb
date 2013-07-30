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

    printf("find_roots....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        fmpcb_poly_t A;
        fmpcb_poly_t B;
        fmpcb_poly_t C;
        fmpcb_t t;
        fmpcb_ptr roots;
        long i, deg, isolated;
        long prec = 10 + n_randint(state, 400);

        fmpcb_init(t);
        fmpcb_poly_init(A);
        fmpcb_poly_init(B);
        fmpcb_poly_init(C);

        do {
            fmpcb_poly_randtest(A, state, 2 + n_randint(state, 15), prec, 5);
        } while (A->length == 0);
        deg = A->length - 1;

        roots = _fmpcb_vec_init(deg);

        isolated = fmpcb_poly_find_roots(roots, A, NULL, 0, prec);

        if (isolated == deg)
        {
            fmpcb_poly_fit_length(B, 1);
            fmpcb_set(B->coeffs, A->coeffs + deg);
            _fmpcb_poly_set_length(B, 1);

            for (i = 0; i < deg; i++)
            {
                fmpcb_poly_fit_length(C, 2);
                fmpcb_one(C->coeffs + 1);
                fmpcb_neg(C->coeffs + 0, roots + i);
                _fmpcb_poly_set_length(C, 2);
                fmpcb_poly_mul(B, B, C, prec);
            }

            if (!fmpcb_poly_contains(B, A))
            {
                printf("FAIL: product does not equal polynomial\n");
                fmpcb_poly_printd(A, 15); printf("\n\n");
                fmpcb_poly_printd(B, 15); printf("\n\n");
                abort();
            }
        }

        for (i = 0; i < isolated; i++)
        {
            fmpcb_poly_evaluate(t, A, roots + i, prec);
            if (!fmpcb_contains_zero(t))
            {
                printf("FAIL: poly(root) does not contain zero\n");
                fmpcb_poly_printd(A, 15); printf("\n\n");
                fmpcb_printd(roots + i, 15); printf("\n\n");
                fmpcb_printd(t, 15); printf("\n\n");
                abort();
            }
        }

        _fmpcb_vec_clear(roots, deg);

        fmpcb_clear(t);
        fmpcb_poly_clear(A);
        fmpcb_poly_clear(B);
        fmpcb_poly_clear(C);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
