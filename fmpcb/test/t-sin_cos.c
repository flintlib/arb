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

    printf("sin_cos....");
    fflush(stdout);

    flint_randinit(state);

    /* check sin(a+b) = cos(b)*sin(a) + cos(a)*sin(b) */
    for (iter = 0; iter < 10000; iter++)
    {
        fmpcb_t a, b, c, d, cosa, sina, cosb, sinb;
        long prec;

        fmpcb_init(a);
        fmpcb_init(b);
        fmpcb_init(c);
        fmpcb_init(d);
        fmpcb_init(cosa);
        fmpcb_init(sina);
        fmpcb_init(cosb);
        fmpcb_init(sinb);

        fmpcb_randtest(a, state, 1 + n_randint(state, 200), 3);
        fmpcb_randtest(b, state, 1 + n_randint(state, 200), 3);

        prec = 2 + n_randint(state, 200);

        fmpcb_add(c, a, b, prec);
        fmpcb_sin(c, c, prec);

        fmpcb_sin_cos(sina, cosa, a, prec);
        fmpcb_sin_cos(sinb, cosb, b, prec);
        fmpcb_mul(cosb, cosb, sina, prec);
        fmpcb_mul(cosa, cosa, sinb, prec);
        fmpcb_add(d, cosa, cosb, prec);

        if (!fmpcb_overlaps(c, d))
        {
            printf("FAIL: sin(a+b) = cos(b)*sin(a) + cos(a)*sin(b)\n\n");
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
        fmpcb_clear(cosa);
        fmpcb_clear(sina);
        fmpcb_clear(cosb);
        fmpcb_clear(sinb);
    }

    /* check cos(a+b) = cos(b)*cos(a) - sin(a)*sin(b) */
    for (iter = 0; iter < 10000; iter++)
    {
        fmpcb_t a, b, c, d, cosa, sina, cosb, sinb;
        long prec;

        fmpcb_init(a);
        fmpcb_init(b);
        fmpcb_init(c);
        fmpcb_init(d);
        fmpcb_init(cosa);
        fmpcb_init(sina);
        fmpcb_init(cosb);
        fmpcb_init(sinb);

        fmpcb_randtest(a, state, 1 + n_randint(state, 200), 3);
        fmpcb_randtest(b, state, 1 + n_randint(state, 200), 3);

        prec = 2 + n_randint(state, 200);

        fmpcb_add(c, a, b, prec);
        fmpcb_cos(c, c, prec);

        fmpcb_sin_cos(sina, cosa, a, prec);
        fmpcb_sin_cos(sinb, cosb, b, prec);
        fmpcb_mul(cosa, cosa, cosb, prec);
        fmpcb_mul(sina, sina, sinb, prec);
        fmpcb_sub(d, cosa, sina, prec);

        if (!fmpcb_overlaps(c, d))
        {
            printf("FAIL: cos(a+b) = cos(b)*cos(a) - sin(a)*sin(b)\n\n");
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
        fmpcb_clear(cosa);
        fmpcb_clear(sina);
        fmpcb_clear(cosb);
        fmpcb_clear(sinb);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
