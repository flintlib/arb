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

    printf("lambda....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        acb_t tau1, tau2, z1, z2, z3, t;
        long e0, prec0, prec1, prec2, step;

        acb_init(tau1);
        acb_init(tau2);
        acb_init(z1);
        acb_init(z2);
        acb_init(z3);
        acb_init(t);

        e0 = 1 + n_randint(state, 100);
        prec0 = 2 + n_randint(state, 1000);
        prec1 = 2 + n_randint(state, 1000);
        prec2 = 2 + n_randint(state, 1000);

        acb_randtest(tau1, state, prec0, e0);
        acb_randtest(tau2, state, prec0, e0);
        acb_randtest(z1, state, prec0, e0);
        acb_randtest(z2, state, prec0, e0);

        acb_set(tau2, tau1);

        step = n_randint(state, 10);

        /* Test invariance */
        while (step --> 0)
        {
            if (n_randint(state, 2))
            {
                acb_add_ui(tau2, tau2, 2, prec0);
            }
            else
            {
                acb_mul_si(z1, tau2, -2, prec0);
                acb_add_ui(z1, z1, 1, prec0);
                acb_div(tau2, tau2, z1, prec0);
            }
        }

        acb_modular_lambda(z1, tau1, prec1);
        acb_modular_lambda(z2, tau2, prec2);

        /* Compare with eta */
        acb_mul_2exp_si(tau1, tau1, -1);
        acb_modular_eta(z3, tau1, prec2);
        acb_mul_2exp_si(tau1, tau1, 2);
        acb_modular_eta(t, tau1, prec2);
        acb_mul(t, t, t, prec2);
        acb_mul(z3, z3, t, prec2);
        acb_mul_2exp_si(tau1, tau1, -1);
        acb_modular_eta(t, tau1, prec2);
        acb_pow_ui(t, t, 3, prec2);
        acb_div(z3, z3, t, prec2);
        acb_pow_ui(z3, z3, 8, prec2);
        acb_mul_2exp_si(z3, z3, 4);

        if (!acb_overlaps(z1, z2) || !acb_overlaps(z1, z3))
        {
            printf("FAIL (overlap)\n");
            printf("tau1 = "); acb_printd(tau1, 15); printf("\n\n");
            printf("tau2 = "); acb_printd(tau2, 15); printf("\n\n");
            printf("z1 = "); acb_printd(z1, 15); printf("\n\n");
            printf("z2 = "); acb_printd(z2, 15); printf("\n\n");
            printf("z3 = "); acb_printd(z3, 15); printf("\n\n");
            abort();
        }

        acb_modular_lambda(tau1, tau1, prec2);

        if (!acb_overlaps(z1, tau1))
        {
            printf("FAIL (aliasing)\n");
            printf("tau1 = "); acb_printd(tau1, 15); printf("\n\n");
            printf("tau2 = "); acb_printd(tau2, 15); printf("\n\n");
            printf("z1 = "); acb_printd(z1, 15); printf("\n\n");
            printf("z2 = "); acb_printd(z2, 15); printf("\n\n");
            abort();
        }

        acb_clear(tau1);
        acb_clear(tau2);
        acb_clear(z1);
        acb_clear(z2);
        acb_clear(z3);
        acb_clear(t);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

