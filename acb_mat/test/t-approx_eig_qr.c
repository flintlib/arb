/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("approx_eig_qr....");
    fflush(stdout);

    flint_randinit(state);

    /* Test DFT matrices */
    for (iter = 0; iter < 50 * arb_test_multiplier(); iter++)
    {
        acb_mat_t A;
        acb_ptr E;
        acb_t t;
        mag_t b;
        slong i, n, prec, goal, c0, c1, c2, c3;

        n = n_randint(state, 30);
        goal = 2 + n_randint(state, 100);

        acb_mat_init(A, n, n);
        acb_init(t);
        mag_init(b);
        E = _acb_vec_init(n);

        for (prec = 32; ; prec *= 2)
        {
            acb_mat_dft(A, 0, prec);
            acb_mat_approx_eig_qr(E, NULL, NULL, A, NULL, 0, prec);

            c0 = c1 = c2 = c3 = 0;

            for (i = 0; i < n; i++)
            {
                acb_set_d_d(t, 1.0, 0.0);
                acb_sub(t, t, E + i, prec);
                acb_get_mag(b, t);
                c0 += (mag_cmp_2exp_si(b, -goal) < 0);

                acb_set_d_d(t, -1.0, 0.0);
                acb_sub(t, t, E + i, prec);
                acb_get_mag(b, t);
                c1 += (mag_cmp_2exp_si(b, -goal) < 0);

                acb_set_d_d(t, 0.0, 1.0);
                acb_sub(t, t, E + i, prec);
                acb_get_mag(b, t);
                c2 += (mag_cmp_2exp_si(b, -goal) < 0);

                acb_set_d_d(t, 0.0, -1.0);
                acb_sub(t, t, E + i, prec);
                acb_get_mag(b, t);
                c3 += (mag_cmp_2exp_si(b, -goal) < 0);
            }

            if (n == 0 || (c0 == (n+4)/4 && c1 == (n+2)/4 && c2 == (n-1)/4 && c3 == (n+1)/4))
                break;

            if (prec > 300)
            {
                flint_printf("FAIL (convergence, DFT matrix)\n\n");
                flint_printf("n = %wd\n\n", n);
                for (i = 0; i < n; i++)
                {
                    acb_printn(E + i, 50, 0);
                    flint_printf("\n");
                }
                flint_abort();
            }
        }

        acb_mat_clear(A);
        _acb_vec_clear(E, n);
        acb_clear(t);
        mag_clear(b);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
