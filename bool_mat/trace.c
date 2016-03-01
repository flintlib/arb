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
bool_mat_trace(const bool_mat_t mat)
{
    slong i;

    if (!bool_mat_is_square(mat))
    {
        flint_printf("bool_mat_trace: a square matrix is required!\n");
        abort();
    }

    if (bool_mat_is_empty(mat))
        return 0;

    for (i = 0; i < bool_mat_nrows(mat); i++)
        if (*bool_mat_entry(mat, i, i))
            return 1;

    return 0;
}
