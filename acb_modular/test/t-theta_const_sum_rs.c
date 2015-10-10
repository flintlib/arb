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

#include "acb_modular.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("theta_const_sum_rs....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 5000; iter++)
    {
        acb_t q, t2a, t2b, t3a, t3b, t4a, t4b;
        long prec1, prec2, N;

        acb_init(q);
        acb_init(t2a);
        acb_init(t2b);
        acb_init(t3a);
        acb_init(t3b);
        acb_init(t4a);
        acb_init(t4b);

        prec1 = 2 + n_randint(state, 1000);
        prec2 = 2 + n_randint(state, 1000);
        N = n_randint(state, 300);

        acb_randtest(q, state, prec1, 3);

        acb_randtest(t2a, state, prec1, 3);
        acb_randtest(t2b, state, prec1, 3);
        acb_randtest(t3a, state, prec1, 3);
        acb_randtest(t3b, state, prec1, 3);
        acb_randtest(t4a, state, prec1, 3);
        acb_randtest(t4b, state, prec1, 3);

        acb_modular_theta_const_sum_basecase(t2a, t3a, t4a, q, N, prec1);
        acb_modular_theta_const_sum_rs(t2b, t3b, t4b, q, N, prec2);

        if (!acb_overlaps(t2a, t2b) || !acb_overlaps(t3a, t3b) || !acb_overlaps(t4a, t4b))
        {
            printf("FAIL (overlap)  iter = %ld\n", iter);
            printf("N = %ld\n", N);
            printf("q = "); acb_printd(q, 50); printf("\n\n");
            printf("t2a = "); acb_printd(t2a, 50); printf("\n\n");
            printf("t2b = "); acb_printd(t2b, 50); printf("\n\n");
            printf("t3a = "); acb_printd(t3a, 50); printf("\n\n");
            printf("t3b = "); acb_printd(t3b, 50); printf("\n\n");
            printf("t4a = "); acb_printd(t4a, 50); printf("\n\n");
            printf("t4b = "); acb_printd(t4b, 50); printf("\n\n");
            abort();
        }

        acb_clear(q);
        acb_clear(t2a);
        acb_clear(t2b);
        acb_clear(t3a);
        acb_clear(t3b);
        acb_clear(t4a);
        acb_clear(t4b);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

