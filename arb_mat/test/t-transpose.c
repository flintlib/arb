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

#include "arb_mat.h"

int main()
{
    slong iter;
    flint_rand_t state;

    printf("transpose....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        slong m, n;
        arb_mat_t a, b, c;

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        arb_mat_init(a, m, n);
        arb_mat_init(b, n, m);
        arb_mat_init(c, m, n);

        arb_mat_randtest(a, state, 2 + n_randint(state, 100), 10);
        arb_mat_randtest(b, state, 2 + n_randint(state, 100), 10);
        arb_mat_randtest(c, state, 2 + n_randint(state, 100), 10);

        arb_mat_transpose(b, a);
        arb_mat_transpose(c, b);

        if (!arb_mat_equal(c, a))
        {
            printf("FAIL\n\n");
            printf("m = %ld, n = %ld\n", m, n);
            abort();
        }

        if (arb_mat_nrows(a) == arb_mat_ncols(a))
        {
            arb_mat_transpose(c, a);
            arb_mat_transpose(a, a);

            if (!arb_mat_equal(a, c))
            {
                printf("FAIL (aliasing)\n\n");
                abort();
            }
        }

        arb_mat_clear(a);
        arb_mat_clear(b);
        arb_mat_clear(c);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
