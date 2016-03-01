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

int
bool_mat_is_transitive(const bool_mat_t mat)
{
    slong n, i, j, k;

    if (!bool_mat_is_square(mat))
    {
        flint_printf("bool_mat_is_transitive: a square matrix is required!\n");
        abort();
    }

    if (bool_mat_is_empty(mat))
        return 1;
    
    n = bool_mat_nrows(mat);

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            for (k = 0; k < n; k++)
                if (*bool_mat_entry(mat, i, j) &&
                    *bool_mat_entry(mat, j, k) &&
                    !*bool_mat_entry(mat, i, k))
                {
                    return 0;
                }

    return 1;
}
