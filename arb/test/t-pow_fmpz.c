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

    printf("pow_fmpz....");
    fflush(stdout);

    flint_randinit(state);

    /* check large arguments */
    for (iter = 0; iter < 10000; iter++)
    {
        arb_t a, b, c, d;
        fmpz_t e1, e2, e3;
        long prec1, prec2;

        prec1 = 2 + n_randint(state, 1000);
        prec2 = prec1 + 30;

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arb_init(d);
        fmpz_init(e1);
        fmpz_init(e2);
        fmpz_init(e3);

        arb_randtest_precise(a, state, 1 + n_randint(state, 1000), 200);
        arb_randtest_precise(b, state, 1 + n_randint(state, 1000), 200);
        fmpz_randtest(e1, state, 200);
        fmpz_randtest(e2, state, 200);

        arb_pow_fmpz(b, a, e1, prec1);
        arb_pow_fmpz(c, a, e1, prec2);

        if (!arb_overlaps(b, c))
        {
            printf("FAIL: overlap\n\n");
            printf("a = "); arb_print(a); printf("\n\n");
            printf("b = "); arb_print(b); printf("\n\n");
            printf("c = "); arb_print(c); printf("\n\n");
            printf("e1 = "); fmpz_print(e1); printf("\n\n");
            abort();
        }

        /* check a^(e1+e2) = a^e1*a^e2 */
        arb_pow_fmpz(c, a, e2, prec1);
        arb_mul(d, b, c, prec1);
        fmpz_add(e3, e1, e2);
        arb_pow_fmpz(c, a, e3, prec1);

        if (!arb_overlaps(c, d))
        {
            printf("FAIL: functional equation\n\n");
            printf("a = "); arb_print(a); printf("\n\n");
            printf("b = "); arb_print(b); printf("\n\n");
            printf("c = "); arb_print(c); printf("\n\n");
            printf("d = "); arb_print(d); printf("\n\n");
            printf("e1 = "); fmpz_print(e1); printf("\n\n");
            printf("e2 = "); fmpz_print(e2); printf("\n\n");
            abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        arb_clear(d);
        fmpz_clear(e1);
        fmpz_clear(e2);
        fmpz_clear(e3);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
