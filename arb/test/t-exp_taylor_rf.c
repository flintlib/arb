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

#include "arb.h"
#include "flint/mpn_extras.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("exp_taylor_rf....");
    fflush(stdout);

    flint_randinit(state);
    _flint_rand_init_gmp(state);

    for (iter = 0; iter < 100000; iter++)
    {
        mp_ptr x, y1, y2, t;
        mp_limb_t err1, err2;
        ulong N;
        mp_size_t xn;
        int cmp, result;

        N = n_randint(state, 288 - 1);
        xn = 1 + n_randint(state, 20);

        x = flint_malloc(sizeof(mp_limb_t) * xn);
        y1 = flint_malloc(sizeof(mp_limb_t) * (xn + 1));
        y2 = flint_malloc(sizeof(mp_limb_t) * (xn + 1));
        t = flint_malloc(sizeof(mp_limb_t) * (xn + 1));

        flint_mpn_rrandom(x, state->gmp_state, xn);
        flint_mpn_rrandom(y1, state->gmp_state, xn + 1);
        flint_mpn_rrandom(y2, state->gmp_state, xn + 1);
        x[xn - 1] &= (LIMB_ONES >> 4);

        _arb_exp_taylor_naive(y1, &err1, x, xn, N);
        _arb_exp_taylor_rs(y2, &err2, x, xn, N);

        cmp = mpn_cmp(y1, y2, xn + 1);

        if (cmp == 0)
        {
            result = 1;
        }
        else if (cmp > 0)
        {
            mpn_sub_n(t, y1, y2, xn + 1);
            result = flint_mpn_zero_p(t + 1, xn) && (t[0] <= err2);
        }
        else
        {
            mpn_sub_n(t, y2, y1, xn + 1);
            result = flint_mpn_zero_p(t + 1, xn) && (t[0] <= err2);
        }

        if (!result)
        {
            flint_printf("FAIL\n");
            flint_printf("N = %wd xn = %wd\n", N, xn);
            flint_printf("x =");
            flint_mpn_debug(x, xn);
            flint_printf("y1 =");
            flint_mpn_debug(y1, xn + 1);
            flint_printf("y2 =");
            flint_mpn_debug(y2, xn + 1);
            abort();
        }

        flint_free(x);
        flint_free(y1);
        flint_free(y2);
        flint_free(t);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

