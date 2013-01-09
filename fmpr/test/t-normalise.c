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

#include "fmpr.h"


int main()
{
    long iter;
    flint_rand_t state;

    printf("normalise....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000000; iter++)
    {
        fmpr_t x, y;
        long prec, ret1, ret2;
        fmpr_rnd_t rnd;

        fmpr_init(x);
        fmpr_init(y);

        fmpz_randtest(fmpr_manref(x), state, 2000);
        fmpz_randtest(fmpr_expref(x), state, 200);

        fmpz_set(fmpr_manref(y), fmpr_manref(x));
        fmpz_set(fmpr_expref(y), fmpr_expref(x));

        switch (n_randint(state, 4))
        {
            case 0: rnd = FMPR_RND_DOWN; break;
            case 1: rnd = FMPR_RND_UP; break;
            case 2: rnd = FMPR_RND_FLOOR; break;
            default: rnd = FMPR_RND_CEIL; break;
        }

        prec = 2 + n_randint(state, 2000);

        ret1 = _fmpr_normalise(fmpr_manref(x), fmpr_expref(x), prec, rnd);
        ret2 = _fmpr_normalise_naive(fmpr_manref(y), fmpr_expref(y), prec, rnd);

        if (!fmpr_equal(x, y) || ret1 != ret2)
        {
            printf("FAIL!\n");
            printf("x = "); fmpr_print(x); printf("\n\n");
            printf("y = "); fmpr_print(y); printf("\n\n");
            printf("ret1 = %ld, ret2 = %ld\n", ret1, ret2);
            abort();
        }

        fmpr_clear(x);
        fmpr_clear(y);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

