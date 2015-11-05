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

    printf("jacobi_p....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        acb_t n, a, b, n1, a1, b1, z, res1, res2, res3, s;
        slong prec;

        acb_init(n);
        acb_init(a);
        acb_init(b);
        acb_init(n1);
        acb_init(a1);
        acb_init(b1);
        acb_init(z);
        acb_init(res1);
        acb_init(res2);
        acb_init(res3);
        acb_init(s);

        prec = 2 + n_randint(state, 300);

        if (n_randint(state, 2))
        {
            acb_set_si(n, n_randint(state, 20) - 10);
            acb_set_si(a, n_randint(state, 20) - 10);
            acb_set_si(b, n_randint(state, 20) - 10);
        }
        else
        {
            acb_randtest_param(n, state, 1 + n_randint(state, 400), 10);
            acb_randtest_param(a, state, 1 + n_randint(state, 400), 10);
            acb_randtest_param(b, state, 1 + n_randint(state, 400), 10);
        }

        acb_randtest_param(z, state, 1 + n_randint(state, 400), 10);

        acb_sub_ui(n1, n, 1, prec);
        acb_sub_ui(a1, a, 1, prec);
        acb_sub_ui(b1, b, 1, prec);

        acb_hypgeom_jacobi_p(res1, n, a, b1, z, prec);
        acb_hypgeom_jacobi_p(res2, n, a1, b, z, 2 + n_randint(state, 300));
        acb_hypgeom_jacobi_p(res3, n1, a, b, z, 2 + n_randint(state, 300));

        acb_sub(s, res1, res2, prec);

        if (!acb_overlaps(s, res3))
        {
            printf("FAIL: consistency\n\n");
            printf("iter = %ld\n\n", iter);
            printf("n = "); acb_printd(n, 30); printf("\n\n");
            printf("a = "); acb_printd(a, 30); printf("\n\n");
            printf("b = "); acb_printd(b, 30); printf("\n\n");
            printf("z = "); acb_printd(z, 30); printf("\n\n");
            printf("res1 = "); acb_printd(res1, 30); printf("\n\n");
            printf("res2 = "); acb_printd(res2, 30); printf("\n\n");
            printf("res3 = "); acb_printd(res3, 30); printf("\n\n");
            printf("s = "); acb_printd(s, 30); printf("\n\n");
            abort();
        }

        acb_clear(n);
        acb_clear(a);
        acb_clear(b);
        acb_clear(n1);
        acb_clear(a1);
        acb_clear(b1);
        acb_clear(z);
        acb_clear(res1);
        acb_clear(res2);
        acb_clear(res3);
        acb_clear(s);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

