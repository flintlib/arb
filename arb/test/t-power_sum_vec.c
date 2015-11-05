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

    Copyright (C) 2015 Fredrik Johansson

******************************************************************************/

#include "arb.h"

int main()
{
    slong iter;
    flint_rand_t state;

    printf("power_sum_vec....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        arb_t a, b, s, t;
        arb_ptr res;
        slong aa, bb, k, n, len;
        slong prec;

        len = n_randint(state, 30);
        prec = 2 + n_randint(state, 500);
        aa = n_randint(state, 50) - 50;
        bb = aa + n_randint(state, 50);

        arb_init(a);
        arb_init(b);
        arb_init(s);
        arb_init(t);
        res = _arb_vec_init(len);

        arb_set_si(a, aa);
        arb_set_si(b, bb);
        arb_power_sum_vec(res, a, b, len, prec);

        for (n = 0; n < len; n++)
        {
            arb_zero(s);
            for (k = aa; k < bb; k++)
            {
                arb_set_si(t, k);
                arb_pow_ui(t, t, n, prec);
                arb_add(s, s, t, prec);
            }

            if (!arb_overlaps(res + n, s))
            {
                printf("FAIL: overlap\n\n");
                printf("a = %ld, b = %ld, n = %ld\n\n", aa, bb, n);
                printf("res = "); arb_printd(res + n, 30); printf("\n\n");
                printf("s = "); arb_printd(s, 30); printf("\n\n");
                abort();
            }
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(s);
        arb_clear(t);
        _arb_vec_clear(res, len);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
