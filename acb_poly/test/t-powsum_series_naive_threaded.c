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

    printf("powsum_series_naive_threaded....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 2000; iter++)
    {
        acb_t s, a, q;
        acb_ptr z1, z2;
        long i, n, len, prec;

        acb_init(s);
        acb_init(a);
        acb_init(q);

        if (n_randint(state, 2))
        {
            acb_randtest(s, state, 1 + n_randint(state, 200), 3);
        }
        else
        {
            arb_set_ui(acb_realref(s), 1);
            arb_mul_2exp_si(acb_realref(s), acb_realref(s), -1);
            arb_randtest(acb_imagref(s), state, 1 + n_randint(state, 200), 4);
        }

        if (n_randint(state, 2))
            acb_one(a);
        else
            acb_randtest(a, state, 1 + n_randint(state, 200), 3);

        if (n_randint(state, 2))
            acb_one(q);
        else
            acb_randtest(q, state, 1 + n_randint(state, 200), 3);

        prec = 2 + n_randint(state, 200);
        n = n_randint(state, 100);
        len = 1 + n_randint(state, 10);

        z1 = _acb_vec_init(len);
        z2 = _acb_vec_init(len);

        _acb_poly_powsum_series_naive(z1, s, a, q, n, len, prec);
        flint_set_num_threads(1 + n_randint(state, 3));
        _acb_poly_powsum_series_naive_threaded(z2, s, a, q, n, len, prec);

        for (i = 0; i < len; i++)
        {
            if (!acb_overlaps(z1 + i, z2 + i))
            {
                printf("FAIL: overlap\n\n");
                printf("iter = %ld\n", iter);
                printf("n = %ld, prec = %ld, len = %ld, i = %ld\n\n", n, prec, len, i);
                printf("s = "); acb_printd(s, prec / 3.33); printf("\n\n");
                printf("a = "); acb_printd(a, prec / 3.33); printf("\n\n");
                printf("q = "); acb_printd(q, prec / 3.33); printf("\n\n");
                printf("z1 = "); acb_printd(z1 + i, prec / 3.33); printf("\n\n");
                printf("z2 = "); acb_printd(z2 + i, prec / 3.33); printf("\n\n");
                abort();
            }
        }

        acb_clear(a);
        acb_clear(s);
        acb_clear(q);
        _acb_vec_clear(z1, len);
        _acb_vec_clear(z2, len);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
