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

    printf("j....");
    fflush(stdout);

    flint_randinit(state);

    /* Test SL2Z invariance */
    for (iter = 0; iter < 10000; iter++)
    {
        acb_t tau1, tau2, z1, z2;
        long e0, prec0, prec1, prec2;
        psl2z_t g;

        psl2z_init(g);
        acb_init(tau1);
        acb_init(tau2);
        acb_init(z1);
        acb_init(z2);

        e0 = 1 + n_randint(state, 100);
        prec0 = 2 + n_randint(state, 2000);
        prec1 = 2 + n_randint(state, 2000);
        prec2 = 2 + n_randint(state, 2000);

        acb_randtest(tau1, state, prec0, e0);
        acb_randtest(tau2, state, prec0, e0);
        acb_randtest(z1, state, prec0, e0);
        acb_randtest(z2, state, prec0, e0);

        psl2z_randtest(g, state, 1 + n_randint(state, 200));

        acb_modular_transform(tau2, g, tau1, prec0);

        acb_modular_j(z1, tau1, prec1);
        acb_modular_j(z2, tau2, prec2);

        if (!acb_overlaps(z1, z2))
        {
            printf("FAIL (overlap)\n");
            printf("tau1 = "); acb_print(tau1); printf("\n\n");
            printf("tau2 = "); acb_print(tau2); printf("\n\n");
            printf("z1 = "); acb_print(z1); printf("\n\n");
            printf("z2 = "); acb_print(z2); printf("\n\n");
            abort();
        }

        acb_modular_j(tau1, tau1, prec2);

        if (!acb_overlaps(z1, tau1))
        {
            printf("FAIL (aliasing)\n");
            printf("tau1 = "); acb_print(tau1); printf("\n\n");
            printf("tau2 = "); acb_print(tau2); printf("\n\n");
            printf("z1 = "); acb_print(z1); printf("\n\n");
            printf("z2 = "); acb_print(z2); printf("\n\n");
            abort();
        }

        acb_clear(tau1);
        acb_clear(tau2);
        acb_clear(z1);
        acb_clear(z2);
        psl2z_clear(g);
    }

    /* Test special values */
    for (iter = 0; iter < 100; iter++)
    {
        acb_t tau, z;
        long prec;

        acb_init(tau);
        acb_init(z);

        prec = 2 + n_randint(state, 2000);

        acb_randtest(z, state, prec, 10);

        acb_onei(tau);
        acb_modular_j(z, tau, prec);
        acb_sub_ui(z, z, 1728, prec);

        if (!acb_contains_zero(z))
        {
            printf("FAIL (value 1)\n");
            printf("tau = "); acb_print(tau); printf("\n\n");
            printf("z = "); acb_print(z); printf("\n\n");
            abort();
        }

        acb_set_ui(tau, 2);
        acb_div_ui(tau, tau, 3, prec);
        acb_exp_pi_i(tau, tau, prec);

        acb_modular_j(z, tau, prec);

        if (!acb_contains_zero(z))
        {
            printf("FAIL (value 2)\n");
            printf("tau = "); acb_print(tau); printf("\n\n");
            printf("z = "); acb_print(z); printf("\n\n");
            abort();
        }

        acb_clear(tau);
        acb_clear(z);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

