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

    printf("0f1....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 2000; iter++)
    {
        acb_t a, z, w1, w2;
        long prec1, prec2;
        int regularized, ebits;

        acb_init(a);
        acb_init(z);
        acb_init(w1);
        acb_init(w2);

        prec1 = 2 + n_randint(state, 700);
        prec2 = 2 + n_randint(state, 700);

        if (n_randint(state, 5) == 0)
            ebits = 100;
        else
            ebits = 10;

        acb_randtest_param(a, state, 1 + n_randint(state, 1000), 1 + n_randint(state, ebits));
        acb_randtest_param(z, state, 1 + n_randint(state, 1000), 1 + n_randint(state, ebits));
        acb_randtest(w1, state, 1 + n_randint(state, 1000), 1 + n_randint(state, ebits));
        acb_randtest(w2, state, 1 + n_randint(state, 1000), 1 + n_randint(state, ebits));
        regularized = n_randint(state, 2);

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_0f1_asymp(w1, a, z, regularized, prec1);
                break;
            case 1:
                acb_hypgeom_0f1_direct(w1, a, z, regularized, prec1);
                break;
            default:
                acb_hypgeom_0f1(w1, a, z, regularized, prec1);
        }

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_0f1_asymp(w2, a, z, regularized, prec2);
                break;
            case 1:
                acb_hypgeom_0f1_direct(w2, a, z, regularized, prec2);
                break;
            default:
                acb_hypgeom_0f1(w2, a, z, regularized, prec2);
        }

        if (!acb_overlaps(w1, w2))
        {
            printf("FAIL: consistency\n\n");
            printf("regularized = %d\n\n", regularized);
            printf("a = "); acb_printd(a, 30); printf("\n\n");
            printf("z = "); acb_printd(z, 30); printf("\n\n");
            printf("w1 = "); acb_printd(w1, 30); printf("\n\n");
            printf("w2 = "); acb_printd(w2, 30); printf("\n\n");
            abort();
        }

        acb_clear(a);
        acb_clear(z);
        acb_clear(w1);
        acb_clear(w2);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

