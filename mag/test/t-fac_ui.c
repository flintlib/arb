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

#include "mag.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("fac_ui....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 2000; iter++)
    {
        fmpr_t x, y;
        fmpz_t f;
        mag_t xb;
        ulong n;

        fmpr_init(x);
        fmpr_init(y);
        fmpz_init(f);
        mag_init(xb);

        mag_randtest_special(xb, state, 80);
        n = n_randtest(state) % 2000;

        mag_fac_ui(xb, n);
        fmpz_fac_ui(f, n);
        fmpr_set_fmpz(x, f);
        mag_get_fmpr(y, xb);

        MAG_CHECK_BITS(xb)

        if (!(fmpr_cmpabs(y, x) >= 0))
        {
            printf("FAIL\n\n");
            printf("n = %lu\n\n", n);
            printf("x = "); fmpr_print(x); printf("\n\n");
            printf("y = "); fmpr_print(y); printf("\n\n");
            abort();
        }

        fmpr_clear(x);
        fmpr_clear(y);
        fmpz_clear(f);
        mag_clear(xb);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

