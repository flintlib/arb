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

#include "arb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("chebyshev_u_ui....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        arb_t a, b, c, d, e;
        ulong n;
        long prec;

        n = n_randtest(state);
        prec = 2 + n_randint(state, 300);

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arb_init(d);
        arb_init(e);

        arb_randtest_precise(a, state, 1 + n_randint(state, 300), 5);
        arb_randtest_precise(c, state, 1 + n_randint(state, 300), 5);

        arb_sin_cos(d, b, a, prec);
        arb_chebyshev_u_ui(c, n, b, prec);
        arb_mul(d, c, d, prec);

        if (n == LIMB_ONES)
            arb_mul_2exp_si(e, a, FLINT_BITS);
        else
            arb_mul_ui(e, a, n + 1, prec);

        arb_sin(e, e, prec);

        if (!arb_overlaps(d, e))
        {
            printf("FAIL: sin(a)*U_n(cos(a)) = sin((n+1)a)\n\n");
            printf("n = %lu\n\n", n);
            printf("a = "); arb_printd(a, 15); printf("\n\n");
            printf("b = "); arb_printd(b, 15); printf("\n\n");
            printf("c = "); arb_printd(c, 15); printf("\n\n");
            printf("d = "); arb_printd(d, 15); printf("\n\n");
            printf("e = "); arb_printd(e, 15); printf("\n\n");
            abort();
        }

        arb_chebyshev_u_ui(b, n, b, prec);

        if (!arb_equal(b, c))
        {
            printf("FAIL: aliasing\n\n");
            printf("n = %lu\n\n", n);
            printf("a = "); arb_printd(a, 15); printf("\n\n");
            printf("b = "); arb_printd(b, 15); printf("\n\n");
            printf("c = "); arb_printd(c, 15); printf("\n\n");
            abort();
        }

        arb_randtest(a, state, 1 + n_randint(state, 300), 5);
        arb_randtest(b, state, 1 + n_randint(state, 300), 5);
        arb_randtest(c, state, 1 + n_randint(state, 300), 5);

        arb_chebyshev_u2_ui(b, c, n, a, prec);
        arb_chebyshev_u_ui(d, n, a, prec);

        if (!arb_overlaps(b, d))
        {
            printf("FAIL: U_n\n\n");
            printf("n = %lu\n\n", n);
            printf("a = "); arb_print(a); printf("\n\n");
            printf("b = "); arb_print(b); printf("\n\n");
            printf("c = "); arb_print(c); printf("\n\n");
            printf("b = "); arb_print(b); printf("\n\n");
            abort();
        }

        if (n == 0)
            arb_zero(d);
        else
            arb_chebyshev_u_ui(d, n - 1, a, prec);

        if (!arb_overlaps(c, d))
        {
            printf("FAIL: U_{n-1}\n\n");
            printf("n = %lu\n\n", n);
            printf("a = "); arb_print(a); printf("\n\n");
            printf("b = "); arb_print(b); printf("\n\n");
            printf("c = "); arb_print(c); printf("\n\n");
            printf("b = "); arb_print(b); printf("\n\n");
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

