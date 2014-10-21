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

    printf("eisenstein....");
    fflush(stdout);

    flint_randinit(state);

    /* Test functional equation */
    for (iter = 0; iter < 10000; iter++)
    {
        acb_t tau1, tau2, t;
        acb_ptr r1, r2;
        long e0, prec0, prec1, prec2, len1, len2, i;
        psl2z_t g;

        psl2z_init(g);
        acb_init(tau1);
        acb_init(tau2);
        acb_init(t);

        e0 = 1 + n_randint(state, 200);
        prec0 = 2 + n_randint(state, 2000);
        prec1 = 2 + n_randint(state, 2000);
        prec2 = 2 + n_randint(state, 2000);
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20);

        r1 = _acb_vec_init(len1);
        r2 = _acb_vec_init(len2);

        acb_randtest(tau1, state, prec0, e0);
        acb_randtest(tau2, state, prec0, e0);

        psl2z_randtest(g, state, 1 + n_randint(state, 200));
        acb_modular_transform(tau2, g, tau1, prec0);

        acb_modular_eisenstein(r1, tau1, len1, prec1);
        acb_modular_eisenstein(r2, tau2, len2, prec2);

        for (i = 0; i < FLINT_MIN(len1, len2); i++)
        {
            acb_mul_fmpz(t, tau1, &g->c, prec1);
            acb_add_fmpz(t, t, &g->d, prec1);
            acb_pow_ui(t, t, 2 * i + 4, prec1);
            acb_mul(t, t, r1 + i, prec1);

            if (!acb_overlaps(t, r2 + i))
            {
                printf("FAIL (overlap)\n");
                printf("tau1 = "); acb_printd(tau1, 15); printf("\n\n");
                printf("tau2 = "); acb_printd(tau2, 15); printf("\n\n");
                printf("g = "); psl2z_print(g); printf("\n\n");
                printf("r1 = "); acb_printd(r1 + i, 15); printf("\n\n");
                printf("r2 = "); acb_printd(r2 + i, 15); printf("\n\n");
                printf("t = "); acb_printd(t, 15); printf("\n\n");
                abort();
            }
        }

        acb_clear(tau1);
        acb_clear(tau2);
        acb_clear(t);
        _acb_vec_clear(r1, len1);
        _acb_vec_clear(r2, len2);
        psl2z_clear(g);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

