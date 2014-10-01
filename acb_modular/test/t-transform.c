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

    printf("transform....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        psl2z_t g;
        acb_t z, w1, w2, t;
        long prec;

        psl2z_init(g);
        acb_init(z);
        acb_init(w1);
        acb_init(w2);
        acb_init(t);

        psl2z_randtest(g, state, n_randint(state, 20));
        acb_randtest(z, state, 1 + n_randint(state, 200), 1 + n_randint(state, 20));
        acb_randtest(w1, state, 1 + n_randint(state, 200), 1 + n_randint(state, 20));
        acb_randtest(w2, state, 1 + n_randint(state, 200), 1 + n_randint(state, 20));

        acb_modular_transform(w1, g, z, 2 + n_randint(state, 200));

        prec = 2 + n_randint(state, 200);

        acb_mul_fmpz(t, z, &g->a, prec);
        acb_add_fmpz(t, t, &g->b, prec);

        acb_mul_fmpz(w2, z, &g->c, prec);
        acb_add_fmpz(w2, w2, &g->d, prec);

        acb_div(w2, t, w2, prec);

        if (!acb_overlaps(w1, w2))
        {
            printf("FAIL\n");
            printf("g = "); psl2z_print(g); printf("\n\n");
            printf("z = "); acb_printd(z, 30); printf("\n\n");
            printf("w1 = "); acb_printd(w1, 30); printf("\n\n");
            printf("w2 = "); acb_printd(w2, 30); printf("\n\n");
            abort();
        }

        psl2z_clear(g);
        acb_clear(z);
        acb_clear(w1);
        acb_clear(w2);
        acb_clear(t);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

