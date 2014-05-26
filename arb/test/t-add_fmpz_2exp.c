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

    Copyright (C) 2012-2014 Fredrik Johansson

******************************************************************************/

#include "arb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("add_fmpz_2exp....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        arb_t a, b, c, d;
        fmpz_t x, e;
        long prec;

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arb_init(d);
        fmpz_init(x);
        fmpz_init(e);

        arb_randtest_special(a, state, 1 + n_randint(state, 2000), 100);
        arb_randtest_special(b, state, 1 + n_randint(state, 2000), 100);
        arb_randtest_special(c, state, 1 + n_randint(state, 2000), 100);
        fmpz_randtest(x, state, 1 + n_randint(state, 2000));
        fmpz_randtest(e, state, 1 + n_randint(state, 200));

        prec = 2 + n_randint(state, 2000);

        arb_set_fmpz_2exp(b, x, e);
        arb_add_fmpz_2exp(c, a, x, e, prec);
        arb_add(d, a, b, prec);

        if (!arb_equal(c, d))
        {
            printf("FAIL\n\n");
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
        fmpz_clear(x);
        fmpz_clear(e);
    }

    /* aliasing */
    for (iter = 0; iter < 10000; iter++)
    {
        arb_t a, c;
        fmpz_t x, e;
        long prec;

        arb_init(a);
        arb_init(c);
        fmpz_init(x);
        fmpz_init(e);

        arb_randtest_special(a, state, 1 + n_randint(state, 2000), 100);
        arb_randtest_special(c, state, 1 + n_randint(state, 2000), 100);
        fmpz_randtest(x, state, 1 + n_randint(state, 2000));
        fmpz_randtest(e, state, 1 + n_randint(state, 200));

        prec = 2 + n_randint(state, 2000);

        arb_add_fmpz_2exp(c, a, x, e, prec);
        arb_add_fmpz_2exp(a, a, x, e, prec);

        if (!arb_equal(a, c))
        {
            printf("FAIL (aliasing)\n\n");
            printf("a = "); arb_print(a); printf("\n\n");
            printf("c = "); arb_print(c); printf("\n\n");
            abort();
        }

        arb_clear(a);
        arb_clear(c);
        fmpz_clear(x);
        fmpz_clear(e);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
