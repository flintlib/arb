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

    printf("eta....");
    fflush(stdout);

    flint_randinit(state);

    /* Test functional equation */
    for (iter = 0; iter < 10000; iter++)
    {
        acb_t tau1, tau2, z1, z2, z3, t;
        fmpq_t arg;
        long e0, prec0, prec1, prec2;
        psl2z_t g;

        psl2z_init(g);
        fmpq_init(arg);
        acb_init(tau1);
        acb_init(tau2);
        acb_init(z1);
        acb_init(z2);
        acb_init(z3);
        acb_init(t);

        e0 = 1 + n_randint(state, 200);
        prec0 = 2 + n_randint(state, 2000);
        prec1 = 2 + n_randint(state, 2000);
        prec2 = 2 + n_randint(state, 2000);

        acb_randtest(tau1, state, prec0, e0);
        acb_randtest(tau2, state, prec0, e0);
        acb_randtest(z1, state, prec0, e0);
        acb_randtest(z2, state, prec0, e0);

        psl2z_randtest(g, state, 1 + n_randint(state, 200));
        acb_modular_transform(tau2, g, tau1, prec0);

        acb_modular_eta(z1, tau1, prec1);
        acb_modular_eta(z2, tau2, prec2);

        /* apply transformation */
        fmpq_set_si(arg, acb_modular_epsilon_arg(g), 12);
        arb_sin_cos_pi_fmpq(acb_imagref(t), acb_realref(t), arg, prec1);
        acb_mul(z3, z1, t, prec1);

        acb_mul_fmpz(t, tau1, &g->c, prec1);
        acb_add_fmpz(t, t, &g->d, prec1);
        acb_sqrt(t, t, prec1);
        acb_mul(z3, z3, t, prec1);

        if (!acb_overlaps(z3, z2))
        {
            printf("FAIL (overlap)\n");
            printf("tau1 = "); acb_printd(tau1, 15); printf("\n\n");
            printf("tau2 = "); acb_printd(tau2, 15); printf("\n\n");
            printf("g = "); psl2z_print(g); printf("\n\n");
            printf("z1 = "); acb_printd(z1, 15); printf("\n\n");
            printf("z2 = "); acb_printd(z2, 15); printf("\n\n");
            printf("z3 = "); acb_printd(z3, 15); printf("\n\n");
            abort();
        }

        acb_modular_eta(tau1, tau1, prec2);

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
        acb_clear(z3);
        acb_clear(t);
        psl2z_clear(g);
        fmpq_clear(arg);
    }

    /* Test special values */
    for (iter = 0; iter < 100; iter++)
    {
        acb_t tau, z;
        arb_t t, u;
        long prec;

        acb_init(tau);
        acb_init(z);
        arb_init(t);
        arb_init(u);

        prec = 2 + n_randint(state, 2000);

        acb_randtest(z, state, prec, 10);

        acb_onei(tau);
        acb_modular_eta(z, tau, prec);

        arb_one(t);
        arb_mul_2exp_si(t, t, -2);
        arb_gamma(t, t, prec);
        arb_const_pi(u, prec);
        arb_root(u, u, 4, prec);
        arb_pow_ui(u, u, 3, prec);
        arb_div(t, t, u, prec);
        arb_mul_2exp_si(t, t, -1);

        if (!arb_overlaps(acb_realref(z), t) ||
            !arb_contains_zero(acb_imagref(z)))
        {
            printf("FAIL (value 1)\n");
            printf("tau = "); acb_print(tau); printf("\n\n");
            printf("z = "); acb_print(z); printf("\n\n");
            abort();
        }

        acb_clear(tau);
        acb_clear(z);
        arb_clear(t);
        arb_clear(u);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

