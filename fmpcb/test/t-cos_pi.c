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

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include "fmpcb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("cos_pi....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        fmpcb_t x, y, a, b, c;
        long prec1, prec2;

        prec1 = 2 + n_randint(state, 1000);
        prec2 = prec1 + 30;

        fmpcb_init(x);
        fmpcb_init(y);
        fmpcb_init(a);
        fmpcb_init(b);
        fmpcb_init(c);

        fmpcb_randtest(x, state, 1 + n_randint(state, 1000), 2 + n_randint(state, 100));

        fmpcb_cos_pi(a, x, prec1);
        fmpcb_cos_pi(b, x, prec2);

        /* check consistency */
        if (!fmpcb_overlaps(a, b))
        {
            printf("FAIL: overlap\n\n");
            printf("x = "); fmpcb_print(x); printf("\n\n");
            printf("a = "); fmpcb_print(a); printf("\n\n");
            printf("b = "); fmpcb_print(b); printf("\n\n");
            abort();
        }

        /* compare with cos */
        fmprb_const_pi(fmpcb_realref(c), prec1);
        fmpcb_mul_fmprb(y, x, fmpcb_realref(c), prec1);
        fmpcb_cos(c, y, prec1);

        if (!fmpcb_overlaps(a, c))
        {
            printf("FAIL: functional equation\n\n");
            printf("x = "); fmpcb_print(x); printf("\n\n");
            printf("y = "); fmpcb_print(y); printf("\n\n");
            printf("a = "); fmpcb_print(a); printf("\n\n");
            printf("c = "); fmpcb_print(c); printf("\n\n");
            abort();
        }

        fmpcb_cos_pi(x, x, prec1);

        if (!fmpcb_overlaps(a, x))
        {
            printf("FAIL: aliasing\n\n");
            printf("a = "); fmpcb_print(a); printf("\n\n");
            printf("x = "); fmpcb_print(x); printf("\n\n");
            abort();
        }

        fmpcb_clear(x);
        fmpcb_clear(y);
        fmpcb_clear(a);
        fmpcb_clear(b);
        fmpcb_clear(c);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

