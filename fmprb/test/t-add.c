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

#include "fmprb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("add....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        fmprb_t a, b, c;
        fmpq_t x, y, z;

        fmprb_init(a);
        fmprb_init(b);
        fmprb_init(c);

        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);

        fmprb_randtest(a, state, 1 + n_randint(state, 200), 10);
        fmprb_randtest(b, state, 1 + n_randint(state, 200), 10);
        fmprb_randtest(c, state, 1 + n_randint(state, 200), 10);

        fmprb_get_rand_fmpq(x, state, a, 1 + n_randint(state, 200));
        fmprb_get_rand_fmpq(y, state, b, 1 + n_randint(state, 200));

        fmprb_add(c, a, b, 2 + n_randint(state, 200));
        fmpq_add(z, x, y);

        if (!fmprb_contains_fmpq(c, z))
        {
            printf("FAIL: containment\n\n");
            printf("a = "); fmprb_print(a); printf("\n\n");
            printf("x = "); fmpq_print(x); printf("\n\n");
            printf("b = "); fmprb_print(b); printf("\n\n");
            printf("y = "); fmpq_print(y); printf("\n\n");
            printf("c = "); fmprb_print(c); printf("\n\n");
            printf("z = "); fmpq_print(z); printf("\n\n");
            abort();
        }

        fmprb_clear(a);
        fmprb_clear(b);
        fmprb_clear(c);

        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_clear(z);
    }

    /* aliasing of c and a */
    for (iter = 0; iter < 10000; iter++)
    {
        fmprb_t a, b;
        fmpq_t x, y, z;

        fmprb_init(a);
        fmprb_init(b);

        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);

        fmprb_randtest(a, state, 1 + n_randint(state, 200), 10);
        fmprb_randtest(b, state, 1 + n_randint(state, 200), 10);

        fmprb_get_rand_fmpq(x, state, a, 1 + n_randint(state, 200));
        fmprb_get_rand_fmpq(y, state, b, 1 + n_randint(state, 200));

        fmprb_add(a, a, b, 2 + n_randint(state, 200));
        fmpq_add(z, x, y);

        if (!fmprb_contains_fmpq(a, z))
        {
            printf("FAIL: aliasing (c, a)\n\n");
            printf("a = "); fmprb_print(a); printf("\n\n");
            printf("x = "); fmpq_print(x); printf("\n\n");
            printf("b = "); fmprb_print(b); printf("\n\n");
            printf("y = "); fmpq_print(y); printf("\n\n");
            printf("z = "); fmpq_print(z); printf("\n\n");
            abort();
        }

        fmprb_clear(a);
        fmprb_clear(b);

        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_clear(z);
    }

    /* aliasing of c and b */
    for (iter = 0; iter < 10000; iter++)
    {
        fmprb_t a, b;
        fmpq_t x, y, z;

        fmprb_init(a);
        fmprb_init(b);

        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);

        fmprb_randtest(a, state, 1 + n_randint(state, 200), 10);
        fmprb_randtest(b, state, 1 + n_randint(state, 200), 10);

        fmprb_get_rand_fmpq(x, state, a, 1 + n_randint(state, 200));
        fmprb_get_rand_fmpq(y, state, b, 1 + n_randint(state, 200));

        fmprb_add(b, a, b, 2 + n_randint(state, 200));
        fmpq_add(z, x, y);

        if (!fmprb_contains_fmpq(b, z))
        {
            printf("FAIL: aliasing (c, b)\n\n");
            printf("a = "); fmprb_print(a); printf("\n\n");
            printf("x = "); fmpq_print(x); printf("\n\n");
            printf("b = "); fmprb_print(b); printf("\n\n");
            printf("y = "); fmpq_print(y); printf("\n\n");
            printf("z = "); fmpq_print(z); printf("\n\n");
            abort();
        }

        fmprb_clear(a);
        fmprb_clear(b);

        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_clear(z);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
