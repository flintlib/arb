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

void
bool_mat_mul_entrywise(bool_mat_t C, const bool_mat_t A, const bool_mat_t B)
{
    slong i, j;

    if (bool_mat_nrows(A) != bool_mat_nrows(B) ||
        bool_mat_ncols(A) != bool_mat_ncols(B))
    {
        flint_printf("bool_mat_mul_entrywise: incompatible dimensions\n");
        abort();
    }

    for (i = 0; i < bool_mat_nrows(A); i++)
    {
        for (j = 0; j < bool_mat_ncols(A); j++)
        {
            bool_mat_set_entry(C, i, j, (bool_mat_get_entry(A, i, j) &
                                         bool_mat_get_entry(B, i, j)));
        }
    }
}
