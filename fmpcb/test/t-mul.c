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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "fmpcb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("mul....");
    fflush(stdout);

    flint_randinit(state);

    /* test aliasing of c and a */
    for (iter = 0; iter < 100000; iter++)
    {
        fmpcb_t a, b, c;
        long prec;

        fmpcb_init(a);
        fmpcb_init(b);
        fmpcb_init(c);

        fmpcb_randtest(a, state, 1 + n_randint(state, 200), 10);
        fmpcb_randtest(b, state, 1 + n_randint(state, 200), 10);
        fmpcb_randtest(c, state, 1 + n_randint(state, 200), 10);

        prec = 2 + n_randint(state, 200);

        fmpcb_mul(c, a, b, prec);
        fmpcb_mul(a, a, b, prec);

        if (!fmpcb_equal(a, c))
        {
            printf("FAIL: aliasing c, a\n\n");
            printf("a = "); fmpcb_print(a); printf("\n\n");
            printf("b = "); fmpcb_print(b); printf("\n\n");
            printf("c = "); fmpcb_print(c); printf("\n\n");
            abort();
        }

        fmpcb_clear(a);
        fmpcb_clear(b);
        fmpcb_clear(c);
    }

    /* test aliasing of c and b */
    for (iter = 0; iter < 100000; iter++)
    {
        fmpcb_t a, b, c;
        long prec;

        fmpcb_init(a);
        fmpcb_init(b);
        fmpcb_init(c);

        fmpcb_randtest(a, state, 1 + n_randint(state, 200), 10);
        fmpcb_randtest(b, state, 1 + n_randint(state, 200), 10);
        fmpcb_randtest(c, state, 1 + n_randint(state, 200), 10);

        prec = 2 + n_randint(state, 200);

        fmpcb_mul(c, a, b, prec);
        fmpcb_mul(b, a, b, prec);

        if (!fmpcb_equal(b, c))
        {
            printf("FAIL: aliasing b, a\n\n");
            printf("a = "); fmpcb_print(a); printf("\n\n");
            printf("b = "); fmpcb_print(b); printf("\n\n");
            printf("c = "); fmpcb_print(c); printf("\n\n");
            abort();
        }

        fmpcb_clear(a);
        fmpcb_clear(b);
        fmpcb_clear(c);
    }

    /* test aliasing a, a */
    for (iter = 0; iter < 100000; iter++)
    {
        fmpcb_t a, b, c, d;
        long prec;

        fmpcb_init(a);
        fmpcb_init(b);
        fmpcb_init(c);
        fmpcb_init(d);

        fmpcb_randtest(a, state, 1 + n_randint(state, 200), 10);
        fmpcb_randtest(b, state, 1 + n_randint(state, 200), 10);
        fmpcb_randtest(c, state, 1 + n_randint(state, 200), 10);

        prec = 2 + n_randint(state, 200);

        fmpcb_set(b, a);
        fmpcb_mul(c, a, a, prec);
        fmpcb_mul(d, a, b, prec);

        if (!fmpcb_equal(c, d))
        {
            printf("FAIL: aliasing a, a\n\n");
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
    }

    /* test aliasing a, a, a */
    for (iter = 0; iter < 100000; iter++)
    {
        fmpcb_t a, b, c;
        long prec;

        fmpcb_init(a);
        fmpcb_init(b);
        fmpcb_init(c);

        fmpcb_randtest(a, state, 1 + n_randint(state, 200), 10);
        fmpcb_randtest(b, state, 1 + n_randint(state, 200), 10);
        fmpcb_randtest(c, state, 1 + n_randint(state, 200), 10);

        prec = 2 + n_randint(state, 200);

        fmpcb_set(b, a);
        fmpcb_mul(c, a, b, prec);
        fmpcb_mul(a, a, a, prec);

        if (!fmpcb_equal(a, c))
        {
            printf("FAIL: aliasing a, a, a\n\n");
            printf("a = "); fmpcb_print(a); printf("\n\n");
            printf("b = "); fmpcb_print(b); printf("\n\n");
            printf("c = "); fmpcb_print(c); printf("\n\n");
            abort();
        }

        fmpcb_clear(a);
        fmpcb_clear(b);
        fmpcb_clear(c);
    }

    /* test a*(b+c) = a*b + a*c */
    for (iter = 0; iter < 100000; iter++)
    {
        fmpcb_t a, b, c, d, e, f;

        fmpcb_init(a);
        fmpcb_init(b);
        fmpcb_init(c);
        fmpcb_init(d);
        fmpcb_init(e);
        fmpcb_init(f);

        fmpcb_randtest(a, state, 1 + n_randint(state, 200), 10);
        fmpcb_randtest(b, state, 1 + n_randint(state, 200), 10);
        fmpcb_randtest(c, state, 1 + n_randint(state, 200), 10);

        fmpcb_add(d, b, c, 2 + n_randint(state, 200));
        fmpcb_mul(e, a, d, 2 + n_randint(state, 200));

        fmpcb_mul(d, a, b, 2 + n_randint(state, 200));
        fmpcb_mul(f, a, c, 2 + n_randint(state, 200));
        fmpcb_add(f, d, f, 2 + n_randint(state, 200));

        if (!fmpcb_overlaps(e, f))
        {
            printf("FAIL: a*(b+c) = a*b + a*c\n\n");
            printf("a = "); fmpcb_print(a); printf("\n\n");
            printf("b = "); fmpcb_print(b); printf("\n\n");
            printf("c = "); fmpcb_print(c); printf("\n\n");
            printf("e = "); fmpcb_print(e); printf("\n\n");
            printf("f = "); fmpcb_print(f); printf("\n\n");
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
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
