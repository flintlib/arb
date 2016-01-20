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
#include "flint/long_extras.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("mul_2exp_si....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        fmpr_t x, y, z;
        mag_t xb, yb;
        slong e;

        fmpr_init(x);
        fmpr_init(y);
        fmpr_init(z);

        mag_init(xb);
        mag_init(yb);

        mag_randtest_special(xb, state, 100);
        e = z_randtest(state);
        mag_get_fmpr(x, xb);

        mag_mul_2exp_si(yb, xb, e);

        fmpr_mul_2exp_si(y, x, e);

        mag_get_fmpr(z, yb);

        MAG_CHECK_BITS(yb)

        if (!fmpr_equal(z, y))
        {
            flint_printf("FAIL\n\n");
            flint_printf("x = "); fmpr_printd(x, 15); flint_printf("\n\n");
            flint_printf("y = "); fmpr_printd(y, 15); flint_printf("\n\n");
            flint_printf("z = "); fmpr_printd(z, 15); flint_printf("\n\n");
            abort();
        }

        fmpr_clear(x);
        fmpr_clear(y);
        fmpr_clear(z);

        mag_clear(xb);
        mag_clear(yb);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
