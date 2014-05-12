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

    Copyright (C) 2012, 2013 Fredrik Johansson

******************************************************************************/

#include "arb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("sin_cos_pi....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        arb_t a, b, c, d, e;
        long prec = 2 + n_randint(state, 200);

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arb_init(d);
        arb_init(e);

        arb_randtest(a, state, 1 + n_randint(state, 200), 10);
        arb_randtest(b, state, 1 + n_randint(state, 200), 10);
        arb_randtest(c, state, 1 + n_randint(state, 200), 10);
        arb_randtest(d, state, 1 + n_randint(state, 200), 10);
        arb_randtest(e, state, 1 + n_randint(state, 200), 10);

        arb_const_pi(b, prec);
        arb_mul(b, b, a, prec);
        arb_sin_cos(b, d, b, prec);

        arb_sin_cos_pi(c, e, a, prec);

        if (!arb_overlaps(b, c) || !arb_overlaps(d, e))
        {
            printf("FAIL: overlap\n\n");
            printf("a = "); arb_print(a); printf("\n\n");
            printf("b = "); arb_print(b); printf("\n\n");
            printf("c = "); arb_print(c); printf("\n\n");
            printf("d = "); arb_print(d); printf("\n\n");
            printf("e = "); arb_print(e); printf("\n\n");
            abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        arb_clear(d);
        arb_clear(e);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

