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

    flint_printf("laguerre_l....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 2000; iter++)
    {
        acb_t n, m, n1, m1, z, res1, res2, res3, s;
        slong prec;

        acb_init(n);
        acb_init(m);
        acb_init(n1);
        acb_init(m1);
        acb_init(z);
        acb_init(res1);
        acb_init(res2);
        acb_init(res3);
        acb_init(s);

        prec = 2 + n_randint(state, 200);

        if (n_randint(state, 2))
        {
            acb_set_si(n, n_randint(state, 20) - 10);
            acb_set_si(m, n_randint(state, 20) - 10);
        }
        else
        {
            acb_randtest_param(n, state, 1 + n_randint(state, 400), 10);
            acb_randtest_param(m, state, 1 + n_randint(state, 400), 10);
        }

        acb_randtest_param(z, state, 1 + n_randint(state, 400), 10);

        acb_sub_ui(n1, n, 1, prec);
        acb_sub_ui(m1, m, 1, prec);

        acb_hypgeom_laguerre_l(res1, n, m, z, prec);
        acb_hypgeom_laguerre_l(res2, n1, m, z, 2 + n_randint(state, 200));
        acb_hypgeom_laguerre_l(res3, n, m1, z, 2 + n_randint(state, 200));

        acb_add(s, res2, res3, prec);

        if (acb_is_finite(res1) && acb_is_finite(s) && !acb_overlaps(res1, s))
        {
            flint_printf("FAIL: consistency\n\n");
            flint_printf("iter = %wd\n\n", iter);
            flint_printf("n = "); acb_printd(n, 30); flint_printf("\n\n");
            flint_printf("m = "); acb_printd(m, 30); flint_printf("\n\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("res1 = "); acb_printd(res1, 30); flint_printf("\n\n");
            flint_printf("res2 = "); acb_printd(res2, 30); flint_printf("\n\n");
            flint_printf("res3 = "); acb_printd(res3, 30); flint_printf("\n\n");
            flint_printf("s = "); acb_printd(s, 30); flint_printf("\n\n");
            abort();
        }

        acb_clear(n);
        acb_clear(m);
        acb_clear(n1);
        acb_clear(m1);
        acb_clear(z);
        acb_clear(res1);
        acb_clear(res2);
        acb_clear(res3);
        acb_clear(s);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

