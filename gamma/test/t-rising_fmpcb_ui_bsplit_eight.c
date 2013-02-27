/*=============================================================================

    This file is part of fmpcb.

    fmpcb is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    fmpcb is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with fmpcb; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "gamma.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("rising_fmpcb_ui_bsplit_eight....");
    fflush(stdout);

    flint_randinit(state);

    /* check functional equation */
    for (iter = 0; iter < 1000; iter++)
    {
        fmpcb_t x, xn, y, z;
        ulong n, m;

        fmpcb_init(x);
        fmpcb_init(xn);
        fmpcb_init(y);
        fmpcb_init(z);

        fmpcb_randtest(x, state, 1 + n_randint(state, 4000), 10);
        n = n_randint(state, 40);
        m = n_randint(state, 40);
        fmpcb_add_ui(xn, x, n, 1 + n_randint(state, 4000));

        gamma_rising_fmpcb_ui_bsplit_eight(y, x, n, 2 + n_randint(state, 4000));
        gamma_rising_fmpcb_ui_bsplit_eight(z, xn, m, 2 + n_randint(state, 4000));
        fmpcb_mul(y, y, z, 2 + n_randint(state, 4000));

        gamma_rising_fmpcb_ui_bsplit_eight(z, x, n + m, 2 + n_randint(state, 4000));

        if (!fmpcb_overlaps(y, z))
        {
            printf("FAIL: overlap\n\n");
            printf("n = %lu\n", n);
            printf("m = %lu\n", m);
            printf("x = "); fmpcb_print(x); printf("\n\n");
            printf("xn = "); fmpcb_print(xn); printf("\n\n");
            printf("y = "); fmpcb_print(y); printf("\n\n");
            printf("z = "); fmpcb_print(z); printf("\n\n");
            abort();
        }

        fmpcb_clear(x);
        fmpcb_clear(xn);
        fmpcb_clear(y);
        fmpcb_clear(z);
    }

    /* aliasing of y and x */
    for (iter = 0; iter < 1000; iter++)
    {
        fmpcb_t x, y;
        ulong n;
        long prec;

        fmpcb_init(x);
        fmpcb_init(y);

        fmpcb_randtest(x, state, 1 + n_randint(state, 200), 10);
        fmpcb_randtest(y, state, 1 + n_randint(state, 200), 10);
        n = n_randint(state, 100);

        prec = 2 + n_randint(state, 4000);
        gamma_rising_fmpcb_ui_bsplit_eight(y, x, n, prec);
        gamma_rising_fmpcb_ui_bsplit_eight(x, x, n, prec);

        if (!fmpcb_equal(x, y))
        {
            printf("FAIL: aliasing\n\n");
            printf("x = "); fmpcb_print(x); printf("\n\n");
            printf("y = "); fmpcb_print(y); printf("\n\n");
            printf("n = %lu\n", n);
            abort();
        }

        fmpcb_clear(x);
        fmpcb_clear(y);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
