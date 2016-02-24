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

#include "fmpz_mat_extras.h"

/* Warshall's algorithm */
void
fmpz_mat_transitive_closure(fmpz_mat_t dest, const fmpz_mat_t src)
{
    slong k, i, j, dim;

    if (fmpz_mat_nrows(dest) != fmpz_mat_nrows(src) ||
        fmpz_mat_ncols(dest) != fmpz_mat_ncols(src))
    {
        flint_printf("fmpz_mat_transitive_closure: incompatible dimensions\n");
        abort();
    }

    dim = fmpz_mat_nrows(src);
    if (dim != fmpz_mat_ncols(src))
    {
        flint_printf("fmpz_mat_transitive_closure: a square matrix is required!\n");
        abort();
    }

    fmpz_mat_set(dest, src);

    for (k = 0; k < dim; k++)
    {
        for (i = 0; i < dim; i++)
        {
            for (j = 0; j < dim; j++)
            {
                if (fmpz_is_zero(fmpz_mat_entry(dest, i, j)) &&
                    !fmpz_is_zero(fmpz_mat_entry(dest, i, k)) &&
                    !fmpz_is_zero(fmpz_mat_entry(dest, k, j)))
                {
                    fmpz_one(fmpz_mat_entry(dest, i, j));
                }
            }
        }
    }
}
