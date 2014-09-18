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

    printf("expm1....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        fmpr_t x, y, z, z2;
        mag_t xb, yb;

        fmpr_init(x);
        fmpr_init(y);
        fmpr_init(z);
        fmpr_init(z2);

        mag_init(xb);
        mag_init(yb);

        mag_randtest_special(xb, state, 0);
        mag_randtest_special(yb, state, 0);

        mag_mul_2exp_si(xb, xb, -100 + n_randint(state,120));

        mag_expm1(yb, xb);

        mag_get_fmpr(x, xb);
        mag_get_fmpr(y, yb);

        fmpr_expm1(z, x, MAG_BITS, FMPR_RND_UP);

        if (fmpr_cmpabs_ui(x, 1000) < 0)
        {
            fmpr_mul_ui(z2, z, 1025, MAG_BITS, FMPR_RND_UP);
            fmpr_mul_2exp_si(z2, z2, -10);
        }
        else
        {
            fmpr_mul_2exp_si(z2, z, 2);
        }

        MAG_CHECK_BITS(xb)
        MAG_CHECK_BITS(yb)

        if (!(fmpr_cmpabs(z, y) <= 0 && fmpr_cmpabs(y, z2) <= 0))
        {
            printf("FAIL\n\n");
            printf("x = "); fmpr_printd(x, 15); printf("\n\n");
            printf("y = "); fmpr_printd(y, 15); printf("\n\n");
            printf("z = "); fmpr_printd(z, 15); printf("\n\n");
            printf("z2 = "); fmpr_printd(z2, 15); printf("\n\n");
            abort();
        }

        mag_expm1(xb, xb);

        if (!mag_equal(xb, yb))
        {
            printf("FAIL (aliasing)\n\n");
            abort();
        }

        fmpr_clear(x);
        fmpr_clear(y);
        fmpr_clear(z);
        fmpr_clear(z2);

        mag_clear(xb);
        mag_clear(yb);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

