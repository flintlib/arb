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

#include "acb_modular.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("psl2z_mul....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        psl2z_t f, g, h, u, v;

        psl2z_init(f);
        psl2z_init(g);
        psl2z_init(h);
        psl2z_init(u);
        psl2z_init(v);

        psl2z_randtest(f, state, n_randint(state, 100));
        psl2z_randtest(g, state, n_randint(state, 100));
        psl2z_randtest(h, state, n_randint(state, 100));
        psl2z_randtest(u, state, n_randint(state, 100));
        psl2z_randtest(v, state, n_randint(state, 100));

        /* test (f*g)*h = f*(g*h) */

        psl2z_mul(u, f, g);
        psl2z_mul(u, u, h);

        psl2z_mul(v, g, h);
        psl2z_mul(v, f, v);

        if (!psl2z_equal(u, v) || !psl2z_is_correct(u) || !psl2z_is_correct(v))
        {
            printf("FAIL\n");
            printf("f = "); psl2z_print(f); printf("\n");
            printf("g = "); psl2z_print(g); printf("\n");
            printf("h = "); psl2z_print(h); printf("\n");
            printf("u = "); psl2z_print(u); printf("\n");
            printf("v = "); psl2z_print(v); printf("\n");
            abort();
        }

        psl2z_clear(f);
        psl2z_clear(g);
        psl2z_clear(h);
        psl2z_clear(u);
        psl2z_clear(v);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

