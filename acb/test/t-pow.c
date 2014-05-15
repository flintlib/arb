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

#include "acb.h"

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
        acb_t a, b, c, d, e, f;
        long prec1, prec2;

        prec1 = 2 + n_randint(state, 1000);
        prec2 = prec1 + 30;

        acb_init(a);
        acb_init(b);
        acb_init(c);
        acb_init(d);
        acb_init(e);
        acb_init(f);

        acb_randtest(a, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 200));
        if (n_randint(state, 4) == 0)
            arb_zero(acb_imagref(a));

        acb_randtest(b, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 200));
        if (n_randint(state, 4) == 0)
            arb_zero(acb_imagref(b));

        acb_pow(c, a, b, prec1);
        acb_pow(d, a, b, prec2);

        if (!acb_overlaps(c, d))
        {
            printf("FAIL: overlap\n\n");
            printf("a = "); acb_print(a); printf("\n\n");
            printf("b = "); acb_print(b); printf("\n\n");
            printf("c = "); acb_print(c); printf("\n\n");
            printf("d = "); acb_print(d); printf("\n\n");
            abort();
        }

        acb_randtest(c, state, 1 + n_randint(state, 1000), 200);

        /* check a^(b+c) = a^b*a^c */
        acb_add(e, b, c, prec1);
        acb_pow(d, a, e, prec1);

        acb_pow(e, a, b, prec1);
        acb_pow(f, a, c, prec1);
        acb_mul(e, e, f, prec1);

        if (!acb_overlaps(d, e))
        {
            printf("FAIL: functional equation\n\n");
            printf("a = "); acb_print(a); printf("\n\n");
            printf("b = "); acb_print(b); printf("\n\n");
            printf("c = "); acb_print(c); printf("\n\n");
            printf("d = "); acb_print(d); printf("\n\n");
            printf("e = "); acb_print(e); printf("\n\n");
            abort();
        }

        acb_pow(c, a, b, prec1);
        acb_set(d, a);
        acb_pow(d, d, b, prec2);

        if (!acb_overlaps(c, d))
        {
            printf("FAIL: aliasing 1\n\n");
            printf("a = "); acb_print(a); printf("\n\n");
            printf("b = "); acb_print(b); printf("\n\n");
            printf("c = "); acb_print(c); printf("\n\n");
            printf("d = "); acb_print(d); printf("\n\n");
            abort();
        }

        acb_set(d, b);
        acb_pow(d, a, d, prec2);

        if (!acb_overlaps(c, d))
        {
            printf("FAIL: aliasing 2\n\n");
            printf("a = "); acb_print(a); printf("\n\n");
            printf("b = "); acb_print(b); printf("\n\n");
            printf("c = "); acb_print(c); printf("\n\n");
            printf("d = "); acb_print(d); printf("\n\n");
            abort();
        }

        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
        acb_clear(d);
        acb_clear(e);
        acb_clear(f);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
