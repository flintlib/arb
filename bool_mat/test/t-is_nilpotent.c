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

    flint_printf("is_nilpotent....");
    fflush(stdout);

    flint_randinit(state);

    /* empty matrix is not nilpotent(?) */
    {
        bool_mat_t A;
        bool_mat_init(A, 0, 0);
        if (bool_mat_is_nilpotent(A))
        {
            flint_printf("FAIL (empty)\n");
            abort();
        }
        bool_mat_clear(A);
    }

    /* detect nilpotency of matrices that are constructed to be nilpotent */
    for (iter = 0; iter < 10000; iter++)
    {
        slong m;
        bool_mat_t A;

        m = n_randint(state, 10) + 1;
        bool_mat_init(A, m, m);
        bool_mat_randtest_nilpotent(A, state);

        if (!bool_mat_is_nilpotent(A))
        {
            flint_printf("FAIL (nilpotent by construction)\n");
            flint_printf("A = \n"); bool_mat_print(A); flint_printf("\n");
            abort();
        }

        bool_mat_clear(A);
    }

    /* check nilpotency by computing a matrix power */
    for (iter = 0; iter < 10000; iter++)
    {
        slong m;
        bool_mat_t A, B;

        m = n_randint(state, 10) + 1;

        bool_mat_init(A, m, m);
        bool_mat_init(B, m, m);

        bool_mat_randtest(A, state);
        bool_mat_pow_ui(B, A, m);

        if (bool_mat_is_nilpotent(A) != !bool_mat_any(B))
        {
            flint_printf("FAIL (A^m)\n");
            flint_printf("A = \n"); bool_mat_print(A); flint_printf("\n");
            flint_printf("B = \n"); bool_mat_print(B); flint_printf("\n");
            abort();
        }

        bool_mat_clear(A);
        bool_mat_clear(B);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
