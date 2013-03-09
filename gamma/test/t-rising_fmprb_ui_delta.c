/*=============================================================================

    This file is part of fmprb.

    fmprb is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    fmprb is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with fmprb; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "gamma.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("rising_fmprb_ui_delta....");
    fflush(stdout);

    flint_randinit(state);

    /* compare with fmpq */
    for (iter = 0; iter < 1000; iter++)
    {
        fmprb_t a, b;
        fmpq_t x, y, z;
        ulong n, step;
        long i;

        fmprb_init(a);
        fmprb_init(b);

        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);

        fmprb_randtest(a, state, 1 + n_randint(state, 1000), 10);
        fmprb_randtest(b, state, 1 + n_randint(state, 1000), 10);
        n = n_randint(state, 80);

        fmprb_get_rand_fmpq(x, state, a, 1 + n_randint(state, 1000));

        step = n_randint(state, 20);
        gamma_rising_fmprb_ui_delta(b, a, n, step, 2 + n_randint(state, 1000));

        fmpq_one(y);
        for (i = 0; i < n; i++)
        {
            fmpq_set_si(z, i, 1);
            fmpq_add(z, x, z);
            fmpq_mul(y, y, z);
        }

        if (!fmprb_contains_fmpq(b, y))
        {
            printf("FAIL: containment\n\n");
            printf("n = %lu\n", n);
            printf("a = "); fmprb_print(a); printf("\n\n");
            printf("x = "); fmpq_print(x); printf("\n\n");
            printf("b = "); fmprb_print(b); printf("\n\n");
            printf("y = "); fmpq_print(y); printf("\n\n");
            abort();
        }

        fmprb_clear(a);
        fmprb_clear(b);

        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_clear(z);
    }

    /* aliasing of y and x */
    for (iter = 0; iter < 1000; iter++)
    {
        fmprb_t x, y;
        ulong n, step;
        long prec;

        fmprb_init(x);
        fmprb_init(y);

        fmprb_randtest(x, state, 1 + n_randint(state, 200), 10);
        fmprb_randtest(y, state, 1 + n_randint(state, 200), 10);
        n = n_randint(state, 100);

        prec = 2 + n_randint(state, 1000);

        step = n_randint(state, 20);
        gamma_rising_fmprb_ui_delta(y, x, n, step, prec);
        gamma_rising_fmprb_ui_delta(x, x, n, step, prec);

        if (!fmprb_equal(x, y))
        {
            printf("FAIL: aliasing\n\n");
            printf("x = "); fmprb_print(x); printf("\n\n");
            printf("y = "); fmprb_print(y); printf("\n\n");
            printf("n = %lu\n", n);
            abort();
        }

        fmprb_clear(x);
        fmprb_clear(y);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
