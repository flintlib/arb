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

    printf("pow....");
    fflush(stdout);

    flint_randinit(state);

    /* check large arguments */
    for (iter = 0; iter < 20000; iter++)
    {
        arb_t a, b, c, d, e, f;
        long prec1, prec2;

        prec1 = 2 + n_randint(state, 1000);
        prec2 = prec1 + 30;

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arb_init(d);
        arb_init(e);
        arb_init(f);

        arb_randtest(a, state, 1 + n_randint(state, 1000), 200);
        arb_randtest(b, state, 1 + n_randint(state, 1000), 200);

        arb_pow(c, a, b, prec1);
        arb_pow(d, a, b, prec2);

        if (!arb_overlaps(c, d))
        {
            printf("FAIL: overlap\n\n");
            printf("a = "); arb_print(a); printf("\n\n");
            printf("b = "); arb_print(b); printf("\n\n");
            printf("c = "); arb_print(c); printf("\n\n");
            printf("d = "); arb_print(d); printf("\n\n");
            abort();
        }

        arb_randtest(c, state, 1 + n_randint(state, 1000), 200);

        /* check a^(b+c) = a^b*a^c */
        arb_add(e, b, c, prec1);
        arb_pow(d, a, e, prec1);

        arb_pow(e, a, b, prec1);
        arb_pow(f, a, c, prec1);
        arb_mul(e, e, f, prec1);

        if (!arb_overlaps(d, e))
        {
            printf("FAIL: functional equation\n\n");
            printf("a = "); arb_print(a); printf("\n\n");
            printf("b = "); arb_print(b); printf("\n\n");
            printf("c = "); arb_print(c); printf("\n\n");
            printf("d = "); arb_print(d); printf("\n\n");
            printf("e = "); arb_print(e); printf("\n\n");
            abort();
        }

        arb_pow(c, a, b, prec1);
        arb_set(d, a);
        arb_pow(d, d, b, prec2);

        if (!arb_overlaps(c, d))
        {
            printf("FAIL: aliasing 1\n\n");
            printf("a = "); arb_print(a); printf("\n\n");
            printf("b = "); arb_print(b); printf("\n\n");
            printf("c = "); arb_print(c); printf("\n\n");
            printf("d = "); arb_print(d); printf("\n\n");
            abort();
        }

        arb_set(d, b);
        arb_pow(d, a, d, prec2);

        if (!arb_overlaps(c, d))
        {
            printf("FAIL: aliasing 2\n\n");
            printf("a = "); arb_print(a); printf("\n\n");
            printf("b = "); arb_print(b); printf("\n\n");
            printf("c = "); arb_print(c); printf("\n\n");
            printf("d = "); arb_print(d); printf("\n\n");
            abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        arb_clear(d);
        arb_clear(e);
        arb_clear(f);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
