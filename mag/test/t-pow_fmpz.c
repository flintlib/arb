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

    printf("pow_fmpz....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        fmpr_t x, y, w;
        mag_t xb, yb;
        fmpz_t e;

        fmpr_init(x);
        fmpr_init(y);
        fmpr_init(w);
        fmpz_init(e);

        mag_init(xb);
        mag_init(yb);

        mag_randtest_special(xb, state, 80);
        mag_randtest_special(yb, state, 80);
        fmpz_randtest(e, state, 200);
        fmpz_abs(e, e);

        mag_get_fmpr(x, xb);

        fmpr_pow_sloppy_fmpz(y, x, e, 2 * MAG_BITS, FMPR_RND_UP);
        if (fmpr_is_nan(y))
            fmpr_pos_inf(y);

        mag_pow_fmpz(yb, xb, e);
        mag_get_fmpr(w, yb);

        MAG_CHECK_BITS(xb)
        MAG_CHECK_BITS(yb)

        if (!(fmpr_cmpabs(y, w) <= 0))
        {
            printf("FAIL\n\n");
            printf("e = "); fmpz_print(e); printf("\n\n");
            printf("x = "); fmpr_print(x); printf("\n\n");
            printf("y = "); fmpr_print(y); printf("\n\n");
            printf("w = "); fmpr_print(w); printf("\n\n");
            abort();
        }

        mag_pow_fmpz(xb, xb, e);

        if (!mag_equal(xb, yb))
        {
            printf("FAIL (aliasing)\n\n");
            printf("e = "); fmpz_print(e); printf("\n\n");
            mag_print(xb); printf("\n\n");
            mag_print(yb); printf("\n\n");
            abort();
        }

        fmpr_clear(x);
        fmpr_clear(y);
        fmpr_clear(w);
        fmpz_clear(e);

        mag_clear(xb);
        mag_clear(yb);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

