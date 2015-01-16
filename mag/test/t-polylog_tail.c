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
#include "arb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("polylog_tail....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 2000; iter++)
    {
        mag_t xb, yb;
        ulong N, k, d;
        long s, prec;
        arb_t z, t, u, r;

        mag_init(xb);
        mag_init(yb);
        arb_init(z);
        arb_init(t);
        arb_init(u);
        arb_init(r);

        mag_randtest_special(xb, state, 6);
        mag_randtest_special(yb, state, 6);
        N = n_randint(state, 100);
        d = n_randint(state, 100);
        s = n_randint(state, 100) - 50;
        prec = 4 * MAG_BITS;

        mag_polylog_tail(yb, xb, s, d, N);

        arb_zero(z);
        arf_set_mag(arb_midref(z), xb);
        arb_zero(r);

        for (k = N; k < N + 100; k++)
        {
            arb_pow_ui(t, z, k, prec);
            arb_log_ui(u, k, prec);
            arb_pow_ui(u, u, d, prec);
            arb_mul(t, t, u, prec);

            arb_set_ui(u, k);
            if (s >= 0)
            {
                arb_pow_ui(u, u, s, prec);
                arb_div(t, t, u, prec);
            }
            else
            {
                arb_pow_ui(u, u, -s, prec);
                arb_mul(t, t, u, prec);
            }

            arb_add(r, r, t, prec);
        }

        MAG_CHECK_BITS(xb)
        MAG_CHECK_BITS(yb)

        arb_zero(z);
        mag_set(arb_radref(z), yb);
        if (!arb_is_finite(z))
            arb_indeterminate(z);

        if (!arb_contains(z, r))
        {
            printf("FAIL\n\n");
            printf("N = %lu\n\n", N);
            printf("d = %lu\n\n", d);
            printf("s = %ld\n\n", s);
            printf("xb = "); mag_printd(xb, 15); printf("\n\n");
            printf("yb = "); mag_printd(yb, 15); printf("\n\n");
            printf("z = "); arb_printd(z, 15); printf("\n\n");
            printf("r = "); arb_printd(r, 15); printf("\n\n");
            abort();
        }

        mag_polylog_tail(xb, xb, s, d, N);

        if (!mag_equal(xb, yb))
        {
            printf("FAIL (aliasing)\n\n");
            abort();
        }

        mag_clear(xb);
        mag_clear(yb);
        arb_clear(z);
        arb_clear(t);
        arb_clear(u);
        arb_clear(r);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

