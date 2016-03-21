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

#include "flint/double_extras.h"
#include "mag.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("get_d....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        mag_t x;
        fmpr_t a, b;
        double y;

        mag_init(x);
        fmpr_init(a);
        fmpr_init(b);

        mag_randtest_special(x, state, 1 + n_randint(state, 100));

        y = mag_get_d(x);
        fmpr_set_d(a, y);

        mag_get_fmpr(b, x);

        if (mag_cmp_2exp_si(x, 1000) < 0 && mag_cmp_2exp_si(x, -1000) > 0)
        {
            if (!fmpr_equal(a, b))
            {
                flint_printf("FAIL (equality)\n\n");
                flint_printf("x = "); mag_print(x); flint_printf("\n\n");
                flint_printf("a = "); fmpr_print(a); flint_printf("\n\n");
                flint_printf("b = "); fmpr_print(b); flint_printf("\n\n");
                abort();
            }
        }
        else
        {
            if (fmpr_cmp(a, b) < 0)
            {
                flint_printf("FAIL (bound)\n\n");
                flint_printf("x = "); mag_print(x); flint_printf("\n\n");
                flint_printf("a = "); fmpr_print(a); flint_printf("\n\n");
                flint_printf("b = "); fmpr_print(b); flint_printf("\n\n");
                abort();
            }
        }

        fmpr_clear(a);
        fmpr_clear(b);
        mag_clear(x);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

