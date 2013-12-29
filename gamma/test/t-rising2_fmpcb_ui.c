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

#include "gamma.h"
#include "arith.h"
#include "fmpcb_poly.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("rising2_fmpcb_ui....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        fmpcb_t a, u, v, u2, v2;
        fmpz *f;
        fmpcb_ptr g;
        ulong n;
        long i, prec;

        fmpcb_init(a);
        fmpcb_init(u);
        fmpcb_init(v);
        fmpcb_init(u2);
        fmpcb_init(v2);

        fmpcb_randtest(a, state, 1 + n_randint(state, 4000), 10);
        fmpcb_randtest(u, state, 1 + n_randint(state, 4000), 10);
        fmpcb_randtest(v, state, 1 + n_randint(state, 4000), 10);
        n = n_randint(state, 120);

        f = _fmpz_vec_init(n + 1);
        g = _fmpcb_vec_init(n + 1);

        prec = 2 + n_randint(state, 4000);
        gamma_rising2_fmpcb_ui(u, v, a, n, prec);

        arith_stirling_number_1u_vec(f, n, n + 1);
        for (i = 0; i <= n; i++)
            fmpcb_set_fmpz(g + i, f + i);
        _fmpcb_poly_evaluate(u2, g, n + 1, a, prec);

        _fmpcb_poly_derivative(g, g, n + 1, prec);
        _fmpcb_poly_evaluate(v2, g, n, a, prec);

        if (!fmpcb_overlaps(u, u2) || !fmpcb_overlaps(v, v2))
        {
            printf("FAIL: overlap\n\n");
            printf("n = %lu\n", n);
            printf("a = "); fmpcb_printd(a, 15); printf("\n\n");
            printf("u = "); fmpcb_printd(u, 15); printf("\n\n");
            printf("u2 = "); fmpcb_printd(u2, 15); printf("\n\n");
            printf("v = "); fmpcb_printd(v, 15); printf("\n\n");
            printf("v2 = "); fmpcb_printd(v2, 15); printf("\n\n");
            abort();
        }

        fmpcb_set(u2, a);
        gamma_rising2_fmpcb_ui(u2, v, u2, n, prec);

        if (!fmpcb_equal(u2, u))
        {
            printf("FAIL: aliasing 1\n\n");
            printf("a = "); fmpcb_printd(a, 15); printf("\n\n");
            printf("u = "); fmpcb_printd(u, 15); printf("\n\n");
            printf("u2 = "); fmpcb_printd(u2, 15); printf("\n\n");
            printf("n = %lu\n", n);
            abort();
        }

        fmpcb_set(v2, a);
        gamma_rising2_fmpcb_ui(u, v2, v2, n, prec);

        if (!fmpcb_equal(v2, v))
        {
            printf("FAIL: aliasing 2\n\n");
            printf("a = "); fmpcb_printd(a, 15); printf("\n\n");
            printf("v = "); fmpcb_printd(v, 15); printf("\n\n");
            printf("v2 = "); fmpcb_printd(v2, 15); printf("\n\n");
            printf("n = %lu\n", n);
            abort();
        }

        fmpcb_clear(a);
        fmpcb_clear(u);
        fmpcb_clear(v);
        fmpcb_clear(u2);
        fmpcb_clear(v2);
        _fmpz_vec_clear(f, n + 1);
        _fmpcb_vec_clear(g, n + 1);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

