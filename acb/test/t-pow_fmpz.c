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

    printf("pow_fmpz....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 5000; iter++)
    {
        acb_t a, b, c, d;
        fmpz_t e1, e2, e3;
        long prec1, prec2;

        prec1 = 2 + n_randint(state, 1000);
        prec2 = prec1 + 30;

        acb_init(a);
        acb_init(b);
        acb_init(c);
        acb_init(d);
        fmpz_init(e1);
        fmpz_init(e2);
        fmpz_init(e3);

        acb_randtest(a, state, 1 + n_randint(state, 1000), 200);
        acb_randtest(b, state, 1 + n_randint(state, 1000), 200);
        fmpz_randtest(e1, state, 200);
        fmpz_randtest(e2, state, 200);

        acb_pow_fmpz(b, a, e1, prec1);
        acb_pow_fmpz(c, a, e1, prec2);

        if (!acb_overlaps(b, c))
        {
            printf("FAIL: overlap\n\n");
            printf("a = "); acb_print(a); printf("\n\n");
            printf("b = "); acb_print(b); printf("\n\n");
            printf("c = "); acb_print(c); printf("\n\n");
            printf("e1 = "); fmpz_print(e1); printf("\n\n");
            abort();
        }

        /* check a^(e1+e2) = a^e1*a^e2 */
        acb_pow_fmpz(c, a, e2, prec1);
        acb_mul(d, b, c, prec1);
        fmpz_add(e3, e1, e2);
        acb_pow_fmpz(c, a, e3, prec1);

        if (!acb_overlaps(c, d))
        {
            printf("FAIL: functional equation\n\n");
            printf("a = "); acb_print(a); printf("\n\n");
            printf("b = "); acb_print(b); printf("\n\n");
            printf("c = "); acb_print(c); printf("\n\n");
            printf("d = "); acb_print(d); printf("\n\n");
            printf("e1 = "); fmpz_print(e1); printf("\n\n");
            printf("e2 = "); fmpz_print(e2); printf("\n\n");
            abort();
        }

        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
        acb_clear(d);
        fmpz_clear(e1);
        fmpz_clear(e2);
        fmpz_clear(e3);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
