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

#include "acb_mat.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("transpose....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        long m, n;
        acb_mat_t a, b, c;

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        acb_mat_init(a, m, n);
        acb_mat_init(b, n, m);
        acb_mat_init(c, m, n);

        acb_mat_randtest(a, state, 2 + n_randint(state, 100), 10);
        acb_mat_randtest(b, state, 2 + n_randint(state, 100), 10);
        acb_mat_randtest(c, state, 2 + n_randint(state, 100), 10);

        acb_mat_transpose(b, a);
        acb_mat_transpose(c, b);

        if (!acb_mat_equal(c, a))
        {
            printf("FAIL\n\n");
            printf("m = %ld, n = %ld\n", m, n);
            abort();
        }

        if (acb_mat_nrows(a) == acb_mat_ncols(a))
        {
            acb_mat_transpose(c, a);
            acb_mat_transpose(a, a);

            if (!acb_mat_equal(a, c))
            {
                printf("FAIL (aliasing)\n\n");
                abort();
            }
        }

        acb_mat_clear(a);
        acb_mat_clear(b);
        acb_mat_clear(c);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
