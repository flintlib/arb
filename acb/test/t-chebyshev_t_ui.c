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

    Copyright (C) 2014 Fredrik Johansson

******************************************************************************/

#include "acb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("chebyshev_t_ui....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        acb_t a, b, c, d;
        ulong n;
        long prec;

        n = n_randtest(state);
        prec = 2 + n_randint(state, 300);

        acb_init(a);
        acb_init(b);
        acb_init(c);
        acb_init(d);

        acb_randtest(a, state, 1 + n_randint(state, 300), 5);
        acb_randtest(c, state, 1 + n_randint(state, 300), 5);

        acb_cos(b, a, prec);
        acb_chebyshev_t_ui(c, n, b, prec);

        acb_mul_ui(d, a, n, prec);
        acb_cos(d, d, prec);

        if (!acb_overlaps(c, d))
        {
            printf("FAIL: c = T_n(cos(a)) = d = cos(n*a)\n\n");
            printf("n = %lu\n\n", n);
            printf("a = "); acb_print(a); printf("\n\n");
            printf("b = "); acb_print(b); printf("\n\n");
            printf("d = "); acb_print(c); printf("\n\n");
            abort();
        }

        acb_chebyshev_t_ui(b, n, b, prec);

        if (!acb_equal(b, c))
        {
            printf("FAIL: aliasing\n\n");
            printf("n = %lu\n\n", n);
            printf("b = "); acb_print(b); printf("\n\n");
            printf("c = "); acb_print(c); printf("\n\n");
            abort();
        }

        acb_randtest(a, state, 1 + n_randint(state, 300), 5);
        acb_randtest(b, state, 1 + n_randint(state, 300), 5);
        acb_randtest(c, state, 1 + n_randint(state, 300), 5);

        acb_chebyshev_t2_ui(b, c, n, a, prec);
        acb_chebyshev_t_ui(d, n, a, prec);

        if (!acb_overlaps(b, d))
        {
            printf("FAIL: T_n\n\n");
            printf("n = %lu\n\n", n);
            printf("a = "); acb_print(a); printf("\n\n");
            printf("b = "); acb_print(b); printf("\n\n");
            printf("c = "); acb_print(c); printf("\n\n");
            printf("b = "); acb_print(b); printf("\n\n");
            abort();
        }

        if (n == 0)
            acb_set(d, a);
        else
            acb_chebyshev_t_ui(d, n - 1, a, prec);

        if (!acb_overlaps(c, d))
        {
            printf("FAIL: T_{n-1}\n\n");
            printf("n = %lu\n\n", n);
            printf("a = "); acb_print(a); printf("\n\n");
            printf("b = "); acb_print(b); printf("\n\n");
            printf("c = "); acb_print(c); printf("\n\n");
            printf("b = "); acb_print(b); printf("\n\n");
            abort();
        }

        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
        acb_clear(d);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

