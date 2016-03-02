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

    Copyright (C) 2016 Arb authors

******************************************************************************/

#include "bool_mat.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("transpose....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        slong m, n;
        bool_mat_t a, b;

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        bool_mat_init(a, m, n);
        bool_mat_init(b, n, m);

        bool_mat_randtest(a, state);
        bool_mat_randtest(b, state);

        bool_mat_transpose(b, a);

        /* involution */
        {
            bool_mat_t c;
            bool_mat_init(c, m, n);
            bool_mat_randtest(c, state);
            bool_mat_transpose(c, b);
            if (!bool_mat_equal(c, a))
            {
                flint_printf("FAIL (involution)\n");
                flint_printf("m = %wd, n = %wd\n", m, n);
                abort();
            }
            bool_mat_clear(c);
        }

        /* aliasing */
        if (bool_mat_is_square(a))
        {
            bool_mat_transpose(a, a);

            if (!bool_mat_equal(a, b))
            {
                flint_printf("FAIL (aliasing)\n");
                abort();
            }
        }

        bool_mat_clear(a);
        bool_mat_clear(b);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
