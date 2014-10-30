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

    Copyright (C) 2014 Fredrik Johansson

******************************************************************************/

#include "arf.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("floor....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        arf_t x, y;
        int result;

        arf_init(x);
        arf_init(y);

        arf_randtest_special(x, state, 2000, 100);
        arf_randtest_special(y, state, 2000, 100);

        arf_floor(y, x);

        result = 1;

        if (arf_is_int(x) || !arf_is_finite(x))
        {
            result = arf_equal(y, x);
        }
        else if (!arf_is_int(y))
        {
            result = 0;
        }
        else if (arf_cmp(y, x) >= 0)
        {
            result = 0;
        }
        else
        {
            arf_t s, t[3];

            /* check floor(x) - x + 1 > 0 */

            arf_init(s);
            arf_init(t[0]);
            arf_init(t[1]);
            arf_init(t[2]);

            arf_set(t[0], y);
            arf_neg(t[1], x);
            arf_one(t[2]);

            arf_sum(s, (arf_ptr) t, 3, 32, ARF_RND_DOWN);

            result = arf_sgn(s) > 0;

            arf_clear(s);
            arf_clear(t[0]);
            arf_clear(t[1]);
            arf_clear(t[2]);
        }

        if (!result)
        {
            printf("FAIL!\n");
            printf("x = "); arf_print(x); printf("\n\n");
            printf("y = "); arf_print(y); printf("\n\n");
            abort();
        }

        arf_floor(x, x);

        if (!arf_equal(x, y))
        {
            printf("FAIL (aliasing)!\n");
            printf("x = "); arf_print(x); printf("\n\n");
            printf("y = "); arf_print(y); printf("\n\n");
            abort();
        }

        arf_clear(x);
        arf_clear(y);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
