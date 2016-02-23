/*=============================================================================

    This file is part of ARB.

    ARB is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    ARB is distributed in the hope that it will be useful,
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

#include "arb_mat.h"

void
_arb_mat_entrywise_nonzero_round_up(fmpz_mat_t A, const arb_mat_t src)
{
    slong i, j;
    fmpz_mat_zero(A);
    for (i = 0; i < arb_mat_nrows(src); i++)
    {
        for (j = 0; j < arb_mat_ncols(src); j++)
        {
            if (!arb_is_zero(arb_mat_entry(src, i, j)))
            {
                fmpz_one(fmpz_mat_entry(A, i, j));
            }
        }
    }
}

slong
_arb_mat_count_nonzero_round_up(const arb_mat_t src)
{
    slong nnz;
    slong i, j;
    nnz = 0;
    for (i = 0; i < arb_mat_nrows(src); i++)
    {
        for (j = 0; j < arb_mat_ncols(src); j++)
        {
            if (!arb_is_zero(arb_mat_entry(src, i, j)))
            {
                nnz++;
            }
        }
    }
    return nnz;
}
