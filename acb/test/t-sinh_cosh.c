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

#include "acb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("sinh_cosh....");
    fflush(stdout);

    flint_randinit(state);

    /* check sinh(a+b) = cosh(b)*sinh(a) + cosh(a)*sinh(b) */
    for (iter = 0; iter < 10000; iter++)
    {
        acb_t a, b, c, d, cosa, sina, cosb, sinb;
        long prec;

        acb_init(a);
        acb_init(b);
        acb_init(c);
        acb_init(d);
        acb_init(cosa);
        acb_init(sina);
        acb_init(cosb);
        acb_init(sinb);

        acb_randtest(a, state, 1 + n_randint(state, 200), 3);
        acb_randtest(b, state, 1 + n_randint(state, 200), 3);

        prec = 2 + n_randint(state, 200);

        acb_add(c, a, b, prec);
        acb_sinh(c, c, prec);

        acb_sinh_cosh(sina, cosa, a, prec);
        acb_sinh_cosh(sinb, cosb, b, prec);
        acb_mul(cosb, cosb, sina, prec);
        acb_mul(cosa, cosa, sinb, prec);
        acb_add(d, cosa, cosb, prec);

        if (!acb_overlaps(c, d))
        {
            printf("FAIL: sinh(a+b) = cosh(b)*sinh(a) + cosh(a)*sinh(b)\n\n");
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
        acb_clear(cosa);
        acb_clear(sina);
        acb_clear(cosb);
        acb_clear(sinb);
    }

    /* check cosh(a+b) = cosh(b)*cosh(a) + sinh(a)*sinh(b) */
    for (iter = 0; iter < 10000; iter++)
    {
        acb_t a, b, c, d, cosa, sina, cosb, sinb;
        long prec;

        acb_init(a);
        acb_init(b);
        acb_init(c);
        acb_init(d);
        acb_init(cosa);
        acb_init(sina);
        acb_init(cosb);
        acb_init(sinb);

        acb_randtest(a, state, 1 + n_randint(state, 200), 3);
        acb_randtest(b, state, 1 + n_randint(state, 200), 3);

        prec = 2 + n_randint(state, 200);

        acb_add(c, a, b, prec);
        acb_cosh(c, c, prec);

        acb_sinh_cosh(sina, cosa, a, prec);
        acb_sinh_cosh(sinb, cosb, b, prec);
        acb_mul(cosa, cosa, cosb, prec);
        acb_mul(sina, sina, sinb, prec);
        acb_add(d, cosa, sina, prec);

        if (!acb_overlaps(c, d))
        {
            printf("FAIL: cosh(a+b) = cosh(b)*cosh(a) + sinh(a)*sinh(b)\n\n");
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
        acb_clear(cosa);
        acb_clear(sina);
        acb_clear(cosb);
        acb_clear(sinb);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
