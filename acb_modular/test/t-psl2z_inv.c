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

    printf("psl2z_inv....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        psl2z_t f, g, h, i;

        psl2z_init(f);
        psl2z_init(g);
        psl2z_init(h);
        psl2z_init(i);

        psl2z_randtest(f, state, n_randint(state, 100));
        psl2z_randtest(g, state, n_randint(state, 100));
        psl2z_randtest(h, state, n_randint(state, 100));

        psl2z_inv(g, f);
        psl2z_mul(h, f, g);
        psl2z_one(i);

        if (!psl2z_equal(h, i) || !psl2z_is_correct(g))
        {
            printf("FAIL\n");
            printf("f = "); psl2z_print(f); printf("\n");
            printf("g = "); psl2z_print(g); printf("\n");
            printf("h = "); psl2z_print(h); printf("\n");
            printf("i = "); psl2z_print(i); printf("\n");
            abort();
        }

        psl2z_inv(f, f);

        if (!psl2z_equal(f, g) || !psl2z_is_correct(f))
        {
            printf("FAIL (aliasing)\n");
            printf("f = "); psl2z_print(f); printf("\n");
            printf("g = "); psl2z_print(g); printf("\n");
            abort();
        }

        psl2z_clear(f);
        psl2z_clear(g);
        psl2z_clear(h);
        psl2z_clear(i);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

