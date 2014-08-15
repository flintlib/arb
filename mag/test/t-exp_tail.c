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

    printf("exp_tail....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        fmpr_t x, t, y, z;
        fmpz_t f;
        mag_t xb, yb;
        ulong N, k;

        fmpr_init(x);
        fmpr_init(t);
        fmpr_init(y);
        fmpr_init(z);
        fmpz_init(f);
        mag_init(xb);
        mag_init(yb);

        mag_randtest_special(xb, state, 6);
        mag_randtest_special(yb, state, 6);
        N = n_randint(state, 100);

        mag_exp_tail(yb, xb, N);

        mag_get_fmpr(x, xb);
        mag_get_fmpr(y, yb);

        fmpr_pow_sloppy_ui(t, x, N, MAG_BITS, FMPR_RND_DOWN);
        fmpz_fac_ui(f, N);
        fmpr_div_fmpz(t, t, f, MAG_BITS, FMPR_RND_DOWN);
        fmpr_set(z, t);

        for (k = 1; k < 50; k++)
        {
            fmpr_mul(t, t, x, MAG_BITS, FMPR_RND_DOWN);
            fmpr_div_ui(t, t, N + k, MAG_BITS, FMPR_RND_DOWN);
            fmpr_add(z, z, t, MAG_BITS, FMPR_RND_DOWN);
        }

        MAG_CHECK_BITS(xb)
        MAG_CHECK_BITS(yb)

        if (!(fmpr_cmpabs(z, y) <= 0))
        {
            printf("FAIL\n\n");
            printf("N = %lu\n\n", N);
            printf("x = "); fmpr_print(x); printf("\n\n");
            printf("y = "); fmpr_print(y); printf("\n\n");
            printf("z = "); fmpr_print(z); printf("\n\n");
            abort();
        }

        mag_exp_tail(xb, xb, N);

        if (!mag_equal(xb, yb))
        {
            printf("FAIL (aliasing)\n\n");
            abort();
        }

        fmpr_clear(x);
        fmpr_clear(t);
        fmpr_clear(y);
        fmpr_clear(z);
        fmpz_clear(f);
        mag_clear(xb);
        mag_clear(yb);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

