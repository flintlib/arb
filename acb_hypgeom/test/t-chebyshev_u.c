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

    Copyright (C) 2015 Fredrik Johansson

******************************************************************************/

#include "acb_hypgeom.h"

int main()
{
    slong iter;
    flint_rand_t state;

    printf("chebyshev_u....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 2000; iter++)
    {
        acb_t n, z, a, b, c, t, res1, res2;
        slong prec;

        acb_init(n);
        acb_init(z);
        acb_init(a);
        acb_init(b);
        acb_init(c);
        acb_init(t);
        acb_init(res1);
        acb_init(res2);

        prec = 2 + n_randint(state, 300);

        acb_randtest_param(n, state, 1 + n_randint(state, 400), 1 + n_randint(state, 10));
        acb_randtest_param(z, state, 1 + n_randint(state, 400), 1 + n_randint(state, 10));
        acb_randtest_param(res1, state, 1 + n_randint(state, 400), 1 + n_randint(state, 10));
        acb_randtest_param(res2, state, 1 + n_randint(state, 400), 1 + n_randint(state, 10));

        acb_hypgeom_chebyshev_u(res1, n, z, prec);

        acb_neg(a, n);
        acb_add_ui(b, n, 2, prec);
        acb_set_ui(c, 3);
        acb_mul_2exp_si(c, c, -1);
        acb_sub_ui(t, z, 1, prec);
        acb_neg(t, t);
        acb_mul_2exp_si(t, t, -1);
        acb_hypgeom_2f1(res2, a, b, c, t, 0, 2 + n_randint(state, 300));
        acb_add_ui(t, n, 1, prec);
        acb_mul(res2, res2, t, prec);

        if (!acb_overlaps(res1, res2))
        {
            printf("FAIL: consistency\n\n");
            printf("iter = %wd\n\n", iter);
            printf("n = "); acb_printd(n, 30); printf("\n\n");
            printf("z = "); acb_printd(z, 30); printf("\n\n");
            printf("res1 = "); acb_printd(res1, 30); printf("\n\n");
            printf("res2 = "); acb_printd(res2, 30); printf("\n\n");
            abort();
        }

        acb_clear(n);
        acb_clear(z);
        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
        acb_clear(t);
        acb_clear(res1);
        acb_clear(res2);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

