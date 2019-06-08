/*
    Copyright (C) 2019 D.H.J. Polymath

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

static int
_arb_vec_overlaps(arb_srcptr a, arb_srcptr b, slong len)
{
    slong i;
    for (i = 0; i < len; i++)
    {
        if (!arb_overlaps(a + i, b + i))
        {
            return 0;
        }
    }
    return 1;
}

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("platt_multieval....");
    fflush(stdout);
    flint_randinit(state);

    for (iter = 0; iter < 10 * arb_test_multiplier(); iter++)
    {
        slong prec;
        ulong A, B, N, J, K;
        slong sigma, Tbits;
        fmpz_t T;
        arb_t h;
        arb_ptr v1, v2;

        /* better but slower limits are in parentheses below */
        prec = 2 + n_randint(state, 300);
        sigma = 1 + 2*(1 + n_randint(state, 100)); /* (200) */
        J = 1 + n_randint(state, 100); /* (10000) */
        K = 1 + n_randint(state, 20); /* (50) */
        A = 1 + n_randint(state, 10);
        B = 1 + n_randint(state, 10); /* (500) */
        if (n_randint(state, 2))
            A *= 2;
        else
            B *= 2;
        N = A*B;

        fmpz_init(T);
        Tbits = 5 + n_randint(state, 15);
        fmpz_set_ui(T, n_randtest_bits(state, Tbits));

        arb_init(h);
        arb_set_si(h, 1 + n_randint(state, 20000));
        arb_div_si(h, h, 1000, prec);

        v1 = _arb_vec_init(N);
        v2 = _arb_vec_init(N);
        acb_dirichlet_platt_scaled_lambda_vec(v1, T, A, B, prec);
        acb_dirichlet_platt_multieval(v2, T, A, B, h, J, K, sigma, prec);

        if (!_arb_vec_overlaps(v1, v2, N))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("iter = %wd  prec = %wd\n\n", iter, prec);
            flint_printf("sigma = %wd\n\n", sigma);
            flint_printf("A = %wu  B = %wu  J = %wu  K = %wu\n\n", A, B, J, K);
            flint_printf("T = "); fmpz_print(T); flint_printf("\n\n");
            flint_printf("h = "); arb_printn(h, 30, 0); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(h);
        fmpz_clear(T);
        _arb_vec_clear(v1, N);
        _arb_vec_clear(v2, N);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
