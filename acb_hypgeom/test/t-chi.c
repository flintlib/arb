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

    Copyright (C) 2015 Fredrik Johansson

******************************************************************************/

#include "acb_hypgeom.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("chi....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        acb_t z0, z1, w0, w1;
        long prec0, prec1;

        acb_init(z0);
        acb_init(z1);
        acb_init(w0);
        acb_init(w1);

        prec0 = 2 + n_randint(state, 1000);
        prec1 = 2 + n_randint(state, 1000);

        acb_randtest(z0, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 10));
        acb_randtest(w0, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 10));
        acb_randtest(w1, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 10));

        acb_set(z1, z0);
        if (n_randint(state, 2))
        {
            acb_add(z1, z1, w0, prec0);
            acb_sub(z1, z1, w0, prec0);
        }

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_chi_2f3(w0, z0, prec0);
                break;
            case 1:
                acb_hypgeom_chi_asymp(w0, z0, prec0);
                break;
            default:
                acb_hypgeom_chi(w0, z0, prec0);
        }

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_chi_2f3(w1, z1, prec1);
                break;
            case 1:
                acb_hypgeom_chi_asymp(w1, z1, prec1);
                break;
            default:
                acb_hypgeom_chi(w1, z1, prec1);
        }

        if (!acb_overlaps(w0, w1))
        {
            printf("FAIL: consistency\n\n");
            printf("z0 = "); acb_printd(z0, 30); printf("\n\n");
            printf("z1 = "); acb_printd(z1, 30); printf("\n\n");
            printf("w0 = "); acb_printd(w0, 30); printf("\n\n");
            printf("w1 = "); acb_printd(w1, 30); printf("\n\n");
            abort();
        }

        acb_clear(z0);
        acb_clear(z1);
        acb_clear(w0);
        acb_clear(w1);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

