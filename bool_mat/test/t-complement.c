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

    flint_printf("complement....");
    fflush(stdout);

    flint_randinit(state);

    /* zero matrix */
    {
        slong m, n;
        for (m = 1; m < 10; m++)
        {
            for (n = 1; n < 10; n++)
            {
                bool_mat_t zero, C;
                bool_mat_init(zero, m, n);
                bool_mat_init(C, m, n);
                bool_mat_zero(zero);
                bool_mat_complement(C, zero);
                if (bool_mat_any(zero) || bool_mat_all(zero) ||
                    !bool_mat_any(C) || !bool_mat_all(C))
                {
                    flint_printf("FAIL (zero matrix)\n");
                    abort();
                }
                bool_mat_clear(zero);
                bool_mat_clear(C);
            }
        }
    }

    /* identity matrix */
    {
        slong n;
        for (n = 2; n < 10; n++)
        {
            bool_mat_t one, C;
            bool_mat_init(one, n, n);
            bool_mat_init(C, n, n);
            bool_mat_one(one);
            bool_mat_complement(C, one);
            if (!bool_mat_any(one) || bool_mat_all(one) ||
                !bool_mat_any(C) || bool_mat_all(C))
            {
                flint_printf("FAIL (identity matrix)\n");
                abort();
            }
            bool_mat_clear(one);
            bool_mat_clear(C);
        }
    }

    /* random matrices */
    for (iter = 0; iter < 10000; iter++)
    {
        slong m, n;
        bool_mat_t A, C;

        m = n_randint(state, 10) + 1;
        n = n_randint(state, 10) + 1;

        bool_mat_init(A, m, n);
        bool_mat_init(C, m, n);

        bool_mat_randtest(A, state);
        bool_mat_randtest(C, state);
        bool_mat_complement(C, A);

        if ((bool_mat_all(A) && bool_mat_any(C)) ||
            (bool_mat_all(C) && bool_mat_any(A)))
        {
            flint_printf("FAIL\n");
            abort();
        }

        /* involution */
        {
            bool_mat_t CC;
            bool_mat_init(CC, m, n);
            bool_mat_randtest(CC, state);
            bool_mat_complement(CC, C);
            if (!bool_mat_equal(A, CC))
            {
                flint_printf("FAIL (involution)\n");
                abort();
            }
            bool_mat_clear(CC);
        }

        /* aliasing */
        bool_mat_complement(A, A);
        if (!bool_mat_equal(A, C))
        {
            flint_printf("FAIL (aliasing)\n");
            abort();
        }

        bool_mat_clear(A);
        bool_mat_clear(C);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
