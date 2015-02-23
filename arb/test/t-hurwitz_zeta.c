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
#include "acb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("hurwitz_zeta....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        arb_t a, b, c;
        acb_t d, e, f;
        long prec;

        prec = 2 + n_randint(state, 300);

        arb_init(a);
        arb_init(b);
        arb_init(c);
        acb_init(d);
        acb_init(e);
        acb_init(f);

        arb_randtest_precise(a, state, 1 + n_randint(state, 300), 5);
        arb_randtest_precise(b, state, 1 + n_randint(state, 300), 5);
        arb_randtest_precise(c, state, 1 + n_randint(state, 300), 5);

        acb_set_arb(d, a);
        acb_set_arb(e, b);

        arb_hurwitz_zeta(c, a, b, prec);
        acb_hurwitz_zeta(f, d, e, prec);

        if (!arb_overlaps(c, acb_realref(f)) ||
            (arb_is_finite(c) && !arb_contains_zero(acb_imagref(f))))
        {
            printf("FAIL: overlap\n\n");
            printf("a = "); arb_printd(a, 15); printf("\n\n");
            printf("b = "); arb_printd(b, 15); printf("\n\n");
            printf("c = "); arb_printd(c, 15); printf("\n\n");
            printf("d = "); acb_printd(d, 15); printf("\n\n");
            printf("e = "); acb_printd(e, 15); printf("\n\n");
            printf("f = "); acb_printd(f, 15); printf("\n\n");
            abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        acb_clear(d);
        acb_clear(e);
        acb_clear(f);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
