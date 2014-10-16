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

#include "profiler.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("delta....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        acb_t tau, z1, z2;
        long e0, prec0, prec1, prec2;

        acb_init(tau);
        acb_init(z1);
        acb_init(z2);

        e0 = 1 + n_randint(state, 10);
        prec0 = 2 + n_randint(state, 1000);
        prec1 = 2 + n_randint(state, 1000);
        prec2 = 2 + n_randint(state, 1000);

        acb_randtest(tau, state, prec0, e0);
        acb_randtest(z1, state, prec0, e0);
        acb_randtest(z2, state, prec0, e0);

        /* compare with eta */
        acb_modular_delta(z1, tau, prec1);

        acb_modular_eta(z2, tau, prec2);
        acb_pow_ui(z2, z2, 24, prec2);

        if (!acb_overlaps(z1, z2))
        {
            printf("FAIL (overlap)\n");
            printf("tau = "); acb_printd(tau, 15); printf("\n\n");
            printf("z1 = "); acb_printd(z1, 15); printf("\n\n");
            printf("z2 = "); acb_printd(z2, 15); printf("\n\n");
            abort();
        }

        acb_modular_delta(tau, tau, prec1);

        if (!acb_overlaps(z1, tau))
        {
            printf("FAIL (aliasing)\n");
            printf("tau = "); acb_printd(tau, 15); printf("\n\n");
            printf("z1 = "); acb_printd(z1, 15); printf("\n\n");
            printf("z2 = "); acb_printd(z2, 15); printf("\n\n");
            abort();
        }

        acb_clear(tau);
        acb_clear(z1);
        acb_clear(z2);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

