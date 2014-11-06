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

#include "arb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("get_mag_lower_nonnegative....");
    fflush(stdout);
    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        arb_t a;
        mag_t m;
        int result;
        arf_struct t[3];
        arf_t s;

        arb_init(a);
        mag_init(m);
        arf_init(s);

        arb_randtest_special(a, state, 200, 1 + n_randint(state, 100));
        arb_get_mag_lower_nonnegative(m, a);
        MAG_CHECK_BITS(m)

        if (arb_contains_nonpositive(a))
        {
            result = mag_is_zero(m);
        }
        else
        {
            arf_init_set_shallow(t + 0, arb_midref(a));
            arf_init_neg_mag_shallow(t + 1, arb_radref(a));
            arf_init_neg_mag_shallow(t + 2, m);
            arf_sum(s, t, 3, 16, ARF_RND_DOWN);
            result = (arf_sgn(s) >= 0);
        }

        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a = "); arb_print(a); printf("\n\n");
            printf("m = "); mag_print(m); printf("\n\n");
            abort();
        }

        arb_clear(a);
        mag_clear(m);
        arf_clear(s);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

