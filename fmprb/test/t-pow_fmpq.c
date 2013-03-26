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

#include "fmprb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("pow_fmpq....");
    fflush(stdout);

    flint_randinit(state);

    /* check large arguments */
    for (iter = 0; iter < 10000; iter++)
    {
        fmprb_t a, b, c, d;
        fmpq_t e1, e2, e3;
        long prec1, prec2;

        prec1 = 2 + n_randint(state, 1000);
        prec2 = prec1 + 30;

        fmprb_init(a);
        fmprb_init(b);
        fmprb_init(c);
        fmprb_init(d);
        fmpq_init(e1);
        fmpq_init(e2);
        fmpq_init(e3);

        fmprb_randtest_precise(a, state, 1 + n_randint(state, 1000), 200);
        fmprb_randtest_precise(b, state, 1 + n_randint(state, 1000), 200);
        fmpq_randtest(e1, state, 200);
        fmpq_randtest(e2, state, 200);

        fmprb_pow_fmpq(b, a, e1, prec1);
        fmprb_pow_fmpq(c, a, e1, prec2);

        if (!fmprb_overlaps(b, c))
        {
            printf("FAIL: overlap\n\n");
            printf("a = "); fmprb_print(a); printf("\n\n");
            printf("b = "); fmprb_print(b); printf("\n\n");
            printf("c = "); fmprb_print(c); printf("\n\n");
            printf("e1 = "); fmpq_print(e1); printf("\n\n");
            abort();
        }

        /* check a^(e1+e2) = a^e1*a^e2 */
        fmprb_pow_fmpq(c, a, e2, prec1);
        fmprb_mul(d, b, c, prec1);
        fmpq_add(e3, e1, e2);
        fmprb_pow_fmpq(c, a, e3, prec1);

        if (!fmprb_overlaps(c, d))
        {
            printf("FAIL: functional equation\n\n");
            printf("a = "); fmprb_print(a); printf("\n\n");
            printf("b = "); fmprb_print(b); printf("\n\n");
            printf("c = "); fmprb_print(c); printf("\n\n");
            printf("d = "); fmprb_print(d); printf("\n\n");
            printf("e1 = "); fmpq_print(e1); printf("\n\n");
            printf("e2 = "); fmpq_print(e2); printf("\n\n");
            abort();
        }

        fmprb_clear(a);
        fmprb_clear(b);
        fmprb_clear(c);
        fmprb_clear(d);
        fmpq_clear(e1);
        fmpq_clear(e2);
        fmpq_clear(e3);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
