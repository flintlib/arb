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
    slong iter;
    flint_rand_t state;

    printf("2f1....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 3000; iter++)
    {
        acb_t a, b, c, z, w1, w2, t;
        slong prec1, prec2;
        int reg1, reg2, ebits;
        int alg1, alg2;

        acb_init(a);
        acb_init(b);
        acb_init(c);
        acb_init(z);
        acb_init(w1);
        acb_init(w2);
        acb_init(t);

        prec1 = 2 + n_randint(state, 300);
        prec2 = 2 + n_randint(state, 300);

        if (n_randint(state, 20) == 0)
            ebits = 30;
        else
            ebits = 5;

        switch (n_randint(state, 3))
        {
            case 0:
                acb_set_si(a, n_randint(state, 500));
                acb_set_si(b, n_randint(state, 500));
                acb_set_si(c, n_randint(state, 10));
                break;
            case 1:
                acb_set_si(a, n_randint(state, 200) - 100);
                acb_set_si(b, n_randint(state, 200) - 100);
                acb_set_si(c, n_randint(state, 200) - 100);
                break;
            default:
                acb_randtest_param(a, state, 1 + n_randint(state, 400), 1 + n_randint(state, ebits));
                acb_randtest_param(b, state, 1 + n_randint(state, 400), 1 + n_randint(state, ebits));
                acb_randtest_param(c, state, 1 + n_randint(state, 400), 1 + n_randint(state, ebits));
        }

        acb_randtest_param(z, state, 1 + n_randint(state, 400), 1 + n_randint(state, ebits));
        acb_randtest(w1, state, 1 + n_randint(state, 400), 1 + n_randint(state, ebits));
        acb_randtest(w2, state, 1 + n_randint(state, 400), 1 + n_randint(state, ebits));

        reg1 = n_randint(state, 2);
        reg2 = n_randint(state, 2);

        alg1 = n_randint(state, 10);
        alg2 = n_randint(state, 10);

        switch (alg1)
        {
            case 0:
                acb_hypgeom_2f1_direct(w1, a, b, c, z, reg1, prec1);
                break;
            case 1:
            case 2:
            case 3:
            case 4:
            case 5:
                acb_hypgeom_2f1_transform(w1, a, b, c, z, reg1, alg1, prec1);
                break;
            case 6:
                acb_hypgeom_2f1_corner(w1, a, b, c, z, reg1, prec1);
                break;
            default:
                acb_hypgeom_2f1(w1, a, b, c, z, reg1, prec1);
        }

        switch (alg2)
        {
            case 0:
                acb_hypgeom_2f1_direct(w2, a, b, c, z, reg2, prec2);
                break;
            case 1:
            case 2:
            case 3:
            case 4:
            case 5:
                acb_hypgeom_2f1_transform(w2, a, b, c, z, reg2, alg2, prec2);
                break;
            case 6:
                acb_hypgeom_2f1_corner(w2, a, b, c, z, reg2, prec2);
                break;
            default:
                acb_hypgeom_2f1(w2, a, b, c, z, reg2, prec2);
        }

        if (reg1 != reg2)
        {
            acb_rgamma(t, c, prec2);

            if (reg1)
                acb_mul(w2, w2, t, prec2);
            else
                acb_mul(w1, w1, t, prec2);
        }

        if (!acb_overlaps(w1, w2))
        {
            printf("FAIL: consistency\n\n");
            printf("iter = %wd, prec1 = %wd, prec2 = %wd\n\n", iter, prec1, prec2);
            printf("alg1 = %d, alg2 = %d\n\n", alg1, alg2);
            printf("reg1 = %d, reg2 = %d\n\n", reg1, reg2);
            printf("a = "); acb_printd(a, 30); printf("\n\n");
            printf("b = "); acb_printd(b, 30); printf("\n\n");
            printf("c = "); acb_printd(c, 30); printf("\n\n");
            printf("z = "); acb_printd(z, 30); printf("\n\n");
            printf("w1 = "); acb_printd(w1, 30); printf("\n\n");
            printf("w2 = "); acb_printd(w2, 30); printf("\n\n");
            abort();
        }

        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
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

