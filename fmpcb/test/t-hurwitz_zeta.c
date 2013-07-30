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

#include "fmpcb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("hurwitz_zeta....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 3000; iter++)
    {
        fmpcb_t a, b, c, d;
        long prec1, prec2;

        prec1 = 2 + n_randint(state, 500);
        prec2 = prec1 + 30;

        fmpcb_init(a);
        fmpcb_init(b);
        fmpcb_init(c);
        fmpcb_init(d);

        fmprb_randtest_precise(fmpcb_realref(a), state, 1 + n_randint(state, 500), 3);
        fmprb_randtest_precise(fmpcb_imagref(a), state, 1 + n_randint(state, 500), 3);

        fmprb_randtest_precise(fmpcb_realref(b), state, 1 + n_randint(state, 500), 3);
        fmprb_randtest_precise(fmpcb_imagref(b), state, 1 + n_randint(state, 500), 3);

        fmpcb_hurwitz_zeta(c, a, b, prec1);
        fmpcb_hurwitz_zeta(d, a, b, prec2);

        if (!fmpcb_overlaps(c, d))
        {
            printf("FAIL: overlap\n\n");
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
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
