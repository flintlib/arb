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

    Copyright (C) 2009 William Hart
    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include "acb_poly.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("shift_left/right....");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing of a and b for left shift */
    for (i = 0; i < 1000; i++)
    {
        acb_poly_t a, b;
        long shift = n_randint(state, 100);

        acb_poly_init(a);
        acb_poly_init(b);
        acb_poly_randtest(a, state, n_randint(state, 100), 2 + n_randint(state, 200), 10);

        acb_poly_shift_left(b, a, shift);
        acb_poly_shift_left(a, a, shift);

        result = (acb_poly_equal(a, b));
        if (!result)
        {
            printf("FAIL:\n");
            acb_poly_printd(a, 10), printf("\n\n");
            acb_poly_printd(b, 10), printf("\n\n");
            abort();
        }

        acb_poly_clear(a);
        acb_poly_clear(b);
    }

    /* Check aliasing of a and b for right shift */
    for (i = 0; i < 1000; i++)
    {
        acb_poly_t a, b;
        long shift = n_randint(state, 100);

        acb_poly_init(a);
        acb_poly_init(b);
        acb_poly_randtest(a, state, n_randint(state, 100), 2 + n_randint(state, 200), 10);

        acb_poly_shift_right(b, a, shift);
        acb_poly_shift_right(a, a, shift);

        result = (acb_poly_equal(a, b));
        if (!result)
        {
            printf("FAIL:\n");
            acb_poly_printd(a, 10), printf("\n\n");
            acb_poly_printd(b, 10), printf("\n\n");
            abort();
        }

        acb_poly_clear(a);
        acb_poly_clear(b);
    }

    /* Check shift left then right does nothing */
    for (i = 0; i < 1000; i++)
    {
        acb_poly_t a, b, c;
        long shift = n_randint(state, 100);

        acb_poly_init(a);
        acb_poly_init(b);
        acb_poly_init(c);
        acb_poly_randtest(a, state, n_randint(state, 100), 2 + n_randint(state, 200), 10);

        acb_poly_shift_left(b, a, shift);
        acb_poly_shift_right(c, b, shift);

        result = (acb_poly_equal(c, a));
        if (!result)
        {
            printf("FAIL:\n");
            acb_poly_printd(a, 10), printf("\n\n");
            acb_poly_printd(b, 10), printf("\n\n");
            acb_poly_printd(c, 10), printf("\n\n");
            abort();
        }

        acb_poly_clear(a);
        acb_poly_clear(b);
        acb_poly_clear(c);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return 0;
}

