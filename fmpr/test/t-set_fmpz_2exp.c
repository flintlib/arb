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

    printf("set_fmpz_2exp....");
    fflush(stdout);

    flint_randinit(state);

    /* test exact roundtrip R -> Q -> R */
    for (iter = 0; iter < 100000; iter++)
    {
        long bits;
        fmpr_t x, z;
        fmpz_t y, e;

        bits = 2 + n_randint(state, 200);

        fmpr_init(x);
        fmpr_init(z);
        fmpz_init(y);
        fmpz_init(e);

        fmpr_randtest(x, state, bits, 10);
        fmpr_randtest(z, state, bits, 10);

        fmpr_get_fmpz_2exp(y, e, x);
        fmpr_set_fmpz_2exp(z, y, e);

        if (!fmpr_equal(x, z))
        {
            printf("FAIL\n\n");
            printf("bits: %ld\n", bits);
            printf("x = "); fmpr_print(x); printf("\n\n");
            printf("y = "); fmpz_print(y); printf("\n\n");
            printf("e = "); fmpz_print(e); printf("\n\n");
            printf("z = "); fmpr_print(z); printf("\n\n");
            abort();
        }

        fmpr_clear(x);
        fmpr_clear(z);
        fmpz_clear(y);
        fmpz_clear(e);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
