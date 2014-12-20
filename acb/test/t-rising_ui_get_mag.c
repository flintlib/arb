/*=============================================================================

    This file is part of acb.

    acb is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    acb is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with acb; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "acb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("rising_ui_get_mag....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        acb_t x, y, z;
        mag_t b;
        ulong n;

        acb_init(x);
        acb_init(y);
        acb_init(z);
        mag_init(b);

        acb_randtest(x, state, 1 + n_randint(state, 400), 1 + n_randint(state, 100));
        n = n_randint(state, 80);

        acb_rising_ui(y, x, n, 2 + n_randint(state, 400));
        acb_rising_ui_get_mag(b, x, n);
        acb_zero(z);
        acb_add_error_mag(z, b);

        if (!acb_overlaps(y, z))
        {
            printf("FAIL: overlap\n\n");
            printf("n = %lu\n", n);
            printf("x = "); acb_printd(x, 15); printf("\n\n");
            printf("y = "); acb_printd(y, 15); printf("\n\n");
            printf("z = "); acb_printd(z, 15); printf("\n\n");
            abort();
        }

        acb_clear(x);
        acb_clear(y);
        acb_clear(z);
        mag_clear(b);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
