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

    printf("theta_1234_sum....");
    fflush(stdout);

    flint_randinit(state);

    /* Very weak test, just testing the error bounds and not
       that we compute the right functions */
    for (iter = 0; iter < 10000; iter++)
    {
        acb_t t1a, t1b, t2a, t2b, t3a, t3b, t4a, t4b, w, q;
        int w_is_unit;
        long prec0, e0, prec1, prec2;

        acb_init(t1a);
        acb_init(t1b);
        acb_init(t2a);
        acb_init(t2b);
        acb_init(t3a);
        acb_init(t3b);
        acb_init(t4a);
        acb_init(t4b);
        acb_init(w);
        acb_init(q);

        e0 = 1 + n_randint(state, 100);
        prec0 = 2 + n_randint(state, 3000);
        prec1 = 2 + n_randint(state, 3000);
        prec2 = 2 + n_randint(state, 3000);

        if (n_randint(state, 2))
        {
            arb_randtest(acb_realref(q), state, prec0, e0);
            arb_zero(acb_imagref(q));
            acb_exp_pi_i(w, q, prec0);
            w_is_unit = n_randint(state, 2);
        }
        else
        {
            acb_randtest(w, state, prec0, e0);
            w_is_unit = 0;
        }

        acb_randtest(q, state, prec0, e0);

        acb_randtest(t1a, state, prec0, e0);
        acb_randtest(t1b, state, prec0, e0);
        acb_randtest(t2a, state, prec0, e0);
        acb_randtest(t2b, state, prec0, e0);
        acb_randtest(t3a, state, prec0, e0);
        acb_randtest(t3b, state, prec0, e0);
        acb_randtest(t4a, state, prec0, e0);
        acb_randtest(t4b, state, prec0, e0);

        acb_modular_theta_1234_sum(t1a, t2a, t3a, t4a, w, w_is_unit, q, 1, prec1);
        acb_modular_theta_1234_sum(t1b, t2b, t3b, t4b, w, w_is_unit & n_randint(state, 2), q, 1, prec2);

        if (!acb_overlaps(t1a, t1b) || !acb_overlaps(t2a, t2b)
            || !acb_overlaps(t3a, t3b) || !acb_overlaps(t4a, t4b))
        {
            printf("FAIL (overlap)\n");
            printf("q = "); acb_print(q); printf("\n\n");
            printf("w = "); acb_print(w); printf("\n\n");
            printf("t1a = "); acb_print(t1a); printf("\n\n");
            printf("t1b = "); acb_print(t1b); printf("\n\n");
            printf("t2a = "); acb_print(t2a); printf("\n\n");
            printf("t2b = "); acb_print(t2b); printf("\n\n");
            printf("t3a = "); acb_print(t3a); printf("\n\n");
            printf("t3b = "); acb_print(t3b); printf("\n\n");
            printf("t4a = "); acb_print(t4a); printf("\n\n");
            printf("t4b = "); acb_print(t4b); printf("\n\n");
            abort();
        }

        acb_clear(t1a);
        acb_clear(t1b);
        acb_clear(t2a);
        acb_clear(t2b);
        acb_clear(t3a);
        acb_clear(t3b);
        acb_clear(t4a);
        acb_clear(t4b);
        acb_clear(w);
        acb_clear(q);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

