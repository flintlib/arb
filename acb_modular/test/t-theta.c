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

    printf("theta....");
    fflush(stdout);

    flint_randinit(state);

    /* Test consistency with/without transform */
    for (iter = 0; iter < 10000; iter++)
    {
        acb_t t1, t2, t3, t4, t1b, t2b, t3b, t4b, z, tau;
        long prec0, prec1, prec2, e0;

        acb_init(t1); acb_init(t2); acb_init(t3); acb_init(t4);
        acb_init(t1b); acb_init(t2b); acb_init(t3b); acb_init(t4b);
        acb_init(z);
        acb_init(tau);

        prec0 = 2 + n_randint(state, 2000);
        prec1 = 2 + n_randint(state, 2000);
        prec2 = 2 + n_randint(state, 2000);
        e0 = 1 + n_randint(state, 100);

        acb_randtest(tau, state, prec0, e0);
        acb_randtest(z, state, prec0, e0);

        acb_modular_theta(t1, t2, t3, t4, z, tau, prec1);
        acb_modular_theta_notransform(t1b, t2b, t3b, t4b, z, tau, prec2);

        if (!acb_overlaps(t1, t1b) || !acb_overlaps(t2, t2b) ||
            !acb_overlaps(t3, t3b) || !acb_overlaps(t4, t4b))
        {
            printf("FAIL (overlap)\n");
            printf("z = "); acb_printd(z, 25); printf("\n\n");
            printf("tau = "); acb_printd(tau, 25); printf("\n\n");

            printf("t1  = "); acb_printd(t1, 15); printf("\n\n");
            printf("t1b = "); acb_printd(t1b, 15); printf("\n\n");

            printf("t2  = "); acb_printd(t2, 15); printf("\n\n");
            printf("t2b = "); acb_printd(t2b, 15); printf("\n\n");

            printf("t3  = "); acb_printd(t3, 15); printf("\n\n");
            printf("t3b = "); acb_printd(t3b, 15); printf("\n\n");

            printf("t4  = "); acb_printd(t4, 15); printf("\n\n");
            printf("t4b = "); acb_printd(t4b, 15); printf("\n\n");

            abort();
        }

        acb_clear(t1); acb_clear(t2); acb_clear(t3); acb_clear(t4);
        acb_clear(t1b); acb_clear(t2b); acb_clear(t3b); acb_clear(t4b);
        acb_clear(z);
        acb_clear(tau);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

