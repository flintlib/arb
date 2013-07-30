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

    printf("exp....");
    fflush(stdout);

    flint_randinit(state);

    /* check exp(a+b) = exp(a)*exp(b) */
    for (iter = 0; iter < 10000; iter++)
    {
        fmpcb_t a, b, c, d, e;
        long prec;

        fmpcb_init(a);
        fmpcb_init(b);
        fmpcb_init(c);
        fmpcb_init(d);
        fmpcb_init(e);

        fmpcb_randtest(a, state, 1 + n_randint(state, 200), 3);
        fmpcb_randtest(b, state, 1 + n_randint(state, 200), 3);
        fmpcb_randtest(c, state, 1 + n_randint(state, 200), 3);

        prec = 2 + n_randint(state, 200);

        fmpcb_add(c, a, b, prec);
        fmpcb_exp(c, c, prec);

        fmpcb_exp(d, a, prec);
        fmpcb_exp(e, b, prec);
        fmpcb_mul(d, d, e, prec);

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
        fmpcb_clear(e);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
