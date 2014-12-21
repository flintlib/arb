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

    printf("binpow_uiui....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        fmpr_t x, y;
        mag_t b;
        ulong m, n;

        fmpr_init(x);
        fmpr_init(y);
        mag_init(b);

        m = n_randtest(state);
        n = n_randtest(state);

        fmpr_one(x);
        fmpr_div_ui(x, x, m, 128, FMPR_RND_UP);
        fmpr_add_ui(x, x, 1, 128, FMPR_RND_UP);
        fmpr_pow_sloppy_ui(x, x, n, 128, FMPR_RND_UP);

        mag_binpow_uiui(b, m, n);
        mag_get_fmpr(y, b);

        MAG_CHECK_BITS(b)

        if (!(fmpr_cmpabs(x, y) <= 0))
        {
            printf("FAIL\n\n");
            printf("m = %lu\n\n", m);
            printf("n = %lu\n\n", n);
            printf("x = "); fmpr_printd(x, 10); printf("\n\n");
            printf("y = "); fmpr_printd(y, 10); printf("\n\n");
            abort();
        }

        fmpr_clear(x);
        fmpr_clear(y);
        mag_clear(b);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

