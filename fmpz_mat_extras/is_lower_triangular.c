/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2016 Arb authors

******************************************************************************/

#include "fmpz_mat_extras.h"

int
fmpz_mat_is_lower_triangular(const fmpz_mat_t mat)
{
    slong i, j;

    if (fmpz_mat_is_empty(mat))
        return 1;

    for (j = 0; j < fmpz_mat_ncols(mat); j++)
    {
        for (i = 0; i < fmpz_mat_nrows(mat) && i < j; i++)
        {
            if (!fmpz_is_zero(fmpz_mat_entry(mat, i, j)))
            {
                return 0;
            }
        }
    }

    return 1;
}
