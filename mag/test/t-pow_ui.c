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
    slong iter;
    flint_rand_t state;

    flint_printf("pow_ui....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        fmpr_t x, y, w;
        mag_t xb, yb;
        ulong e;

        fmpr_init(x);
        fmpr_init(y);
        fmpr_init(w);

        mag_init(xb);
        mag_init(yb);

        mag_randtest_special(xb, state, 80);
        mag_randtest_special(yb, state, 80);
        e = n_randtest(state);

        mag_get_fmpr(x, xb);

        fmpr_pow_sloppy_ui(y, x, e, 2 * MAG_BITS, FMPR_RND_UP);
        if (fmpr_is_nan(y))
            fmpr_pos_inf(y);

        mag_pow_ui(yb, xb, e);
        mag_get_fmpr(w, yb);

        MAG_CHECK_BITS(xb)
        MAG_CHECK_BITS(yb)

        if (!(fmpr_cmpabs(y, w) <= 0))
        {
            flint_printf("FAIL\n\n");
            flint_printf("e = %wu\n\n", e);
            flint_printf("x = "); fmpr_print(x); flint_printf("\n\n");
            flint_printf("y = "); fmpr_print(y); flint_printf("\n\n");
            flint_printf("w = "); fmpr_print(w); flint_printf("\n\n");
            abort();
        }

        mag_pow_ui(xb, xb, e);

        if (!mag_equal(xb, yb))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_printf("e = %wu\n\n", e);
            mag_print(xb); flint_printf("\n\n");
            mag_print(yb); flint_printf("\n\n");
            abort();
        }

        fmpr_clear(x);
        fmpr_clear(y);
        fmpr_clear(w);

        mag_clear(xb);
        mag_clear(yb);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

