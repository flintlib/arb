/*=============================================================================

    This file is part of fmprb.

    fmprb is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    fmprb is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with fmprb; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include "gamma.h"
#include "arith.h"
#include "fmprb_poly.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("rising2_fmprb_ui....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        fmprb_t a, u, v, u2, v2;
        fmpz *f;
        fmprb_ptr g;
        ulong n;
        long i, prec;

        fmprb_init(a);
        fmprb_init(u);
        fmprb_init(v);
        fmprb_init(u2);
        fmprb_init(v2);

        fmprb_randtest(a, state, 1 + n_randint(state, 4000), 10);
        fmprb_randtest(u, state, 1 + n_randint(state, 4000), 10);
        fmprb_randtest(v, state, 1 + n_randint(state, 4000), 10);
        n = n_randint(state, 120);

        f = _fmpz_vec_init(n + 1);
        g = _fmprb_vec_init(n + 1);

        prec = 2 + n_randint(state, 4000);
        gamma_rising2_fmprb_ui(u, v, a, n, prec);

        arith_stirling_number_1u_vec(f, n, n + 1);
        for (i = 0; i <= n; i++)
            fmprb_set_fmpz(g + i, f + i);
        _fmprb_poly_evaluate(u2, g, n + 1, a, prec);

        _fmprb_poly_derivative(g, g, n + 1, prec);
        _fmprb_poly_evaluate(v2, g, n, a, prec);

        if (!fmprb_overlaps(u, u2) || !fmprb_overlaps(v, v2))
        {
            printf("FAIL: overlap\n\n");
            printf("n = %lu\n", n);
            printf("a = "); fmprb_printd(a, 15); printf("\n\n");
            printf("u = "); fmprb_printd(u, 15); printf("\n\n");
            printf("u2 = "); fmprb_printd(u2, 15); printf("\n\n");
            printf("v = "); fmprb_printd(v, 15); printf("\n\n");
            printf("v2 = "); fmprb_printd(v2, 15); printf("\n\n");
            abort();
        }

        fmprb_set(u2, a);
        gamma_rising2_fmprb_ui(u2, v, u2, n, prec);

        if (!fmprb_equal(u2, u))
        {
            printf("FAIL: aliasing 1\n\n");
            printf("a = "); fmprb_printd(a, 15); printf("\n\n");
            printf("u = "); fmprb_printd(u, 15); printf("\n\n");
            printf("u2 = "); fmprb_printd(u2, 15); printf("\n\n");
            printf("n = %lu\n", n);
            abort();
        }

        fmprb_set(v2, a);
        gamma_rising2_fmprb_ui(u, v2, v2, n, prec);

        if (!fmprb_equal(v2, v))
        {
            printf("FAIL: aliasing 2\n\n");
            printf("a = "); fmprb_printd(a, 15); printf("\n\n");
            printf("v = "); fmprb_printd(v, 15); printf("\n\n");
            printf("v2 = "); fmprb_printd(v2, 15); printf("\n\n");
            printf("n = %lu\n", n);
            abort();
        }

        fmprb_clear(a);
        fmprb_clear(u);
        fmprb_clear(v);
        fmprb_clear(u2);
        fmprb_clear(v2);
        _fmpz_vec_clear(f, n + 1);
        _fmprb_vec_clear(g, n + 1);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

