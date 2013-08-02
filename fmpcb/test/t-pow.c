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

#include "fmpcb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("pow....");
    fflush(stdout);

    flint_randinit(state);

    /* check large arguments */
    for (iter = 0; iter < 10000; iter++)
    {
        fmpcb_t a, b, c, d, e, f;
        long prec1, prec2;

        prec1 = 2 + n_randint(state, 1000);
        prec2 = prec1 + 30;

        fmpcb_init(a);
        fmpcb_init(b);
        fmpcb_init(c);
        fmpcb_init(d);
        fmpcb_init(e);
        fmpcb_init(f);

        fmpcb_randtest(a, state, 1 + n_randint(state, 1000), 200);
        fmpcb_randtest(b, state, 1 + n_randint(state, 1000), 200);

        fmpcb_pow(c, a, b, prec1);
        fmpcb_pow(d, a, b, prec2);

        if (!fmpcb_overlaps(c, d))
        {
            printf("FAIL: overlap\n\n");
            printf("a = "); fmpcb_print(a); printf("\n\n");
            printf("b = "); fmpcb_print(b); printf("\n\n");
            printf("c = "); fmpcb_print(c); printf("\n\n");
            printf("d = "); fmpcb_print(d); printf("\n\n");
            abort();
        }

        fmpcb_randtest(c, state, 1 + n_randint(state, 1000), 200);

        /* check a^(b+c) = a^b*a^c */
        fmpcb_add(e, b, c, prec1);
        fmpcb_pow(d, a, e, prec1);

        fmpcb_pow(e, a, b, prec1);
        fmpcb_pow(f, a, c, prec1);
        fmpcb_mul(e, e, f, prec1);

        if (!fmpcb_overlaps(d, e))
        {
            printf("FAIL: functional equation\n\n");
            printf("a = "); fmpcb_print(a); printf("\n\n");
            printf("b = "); fmpcb_print(b); printf("\n\n");
            printf("c = "); fmpcb_print(c); printf("\n\n");
            printf("d = "); fmpcb_print(d); printf("\n\n");
            printf("e = "); fmpcb_print(e); printf("\n\n");
            abort();
        }

        fmpcb_pow(c, a, b, prec1);
        fmpcb_set(d, a);
        fmpcb_pow(d, d, b, prec2);

        if (!fmpcb_overlaps(c, d))
        {
            printf("FAIL: aliasing 1\n\n");
            printf("a = "); fmpcb_print(a); printf("\n\n");
            printf("b = "); fmpcb_print(b); printf("\n\n");
            printf("c = "); fmpcb_print(c); printf("\n\n");
            printf("d = "); fmpcb_print(d); printf("\n\n");
            abort();
        }

        fmpcb_set(d, b);
        fmpcb_pow(d, a, d, prec2);

        if (!fmpcb_overlaps(c, d))
        {
            printf("FAIL: aliasing 2\n\n");
            printf("a = "); fmpcb_print(a); printf("\n\n");
            printf("b = "); fmpcb_print(b); printf("\n\n");
            printf("c = "); fmpcb_print(c); printf("\n\n");
            printf("d = "); fmpcb_print(d); printf("\n\n");
            abort();
        }

        fmpcb_clear(a);
        fmpcb_clear(b);
        fmpcb_clear(c);
        fmpcb_clear(d);
        fmpcb_clear(e);
        fmpcb_clear(f);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
