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

    printf("set_round_uiui....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000000; iter++)
    {
        arf_t x, y;
        long prec, fix1, fix2;
        int ret1, ret2, sgnbit;
        mp_limb_t t[2];
        arf_rnd_t rnd;

        prec = 2 + n_randint(state, 1000);

        arf_init(x);
        arf_init(y);

        arf_randtest_special(x, state, 1 + n_randint(state, 200), 1 + n_randint(state, 100));
        arf_randtest_special(y, state, 1 + n_randint(state, 200), 1 + n_randint(state, 100));

        do {
            t[0] = n_randtest(state);
            t[1] = n_randtest(state);
        } while (t[0] == 0 && t[1] == 0);

        sgnbit = n_randint(state, 2);

        switch (n_randint(state, 4))
        {
            case 0: rnd = ARF_RND_DOWN; break;
            case 1: rnd = ARF_RND_UP; break;
            case 2: rnd = ARF_RND_FLOOR; break;
            default: rnd = ARF_RND_CEIL; break;
        }

        if (t[1] != 0)
        {
            ret1 = _arf_set_round_mpn(x, &fix1, t, 2, sgnbit, prec, rnd);
            fmpz_set_si(ARF_EXPREF(x), 2 * FLINT_BITS + fix1);
        }
        else
        {
            ret1 = _arf_set_round_mpn(x, &fix1, t, 1, sgnbit, prec, rnd);
            fmpz_set_si(ARF_EXPREF(x), FLINT_BITS + fix1);
        }

        ret2 = _arf_set_round_uiui(y, &fix2, t[1], t[0], sgnbit, prec, rnd);
        fmpz_set_si(ARF_EXPREF(y), 2 * FLINT_BITS + fix2);

        if (!arf_equal(x, y) || (ret1 != ret2))
        {
            printf("FAIL\n\n");
            printf("prec = %ld", prec); printf("\n\n");
            printf("hi = %lu, lo = %lu\n\n", t[1], t[0]);
            printf("x = "); arf_print(x); printf("\n\n");
            printf("y = "); arf_print(y); printf("\n\n");
            printf("ret1 = %d, ret2 = %d\n\n", ret1, ret2);
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

