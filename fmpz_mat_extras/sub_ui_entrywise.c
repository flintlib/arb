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

void
fmpz_mat_sub_ui_entrywise(fmpz_mat_t B, const fmpz_mat_t A, ulong x)
{
    slong i, j;
    for (i = 0; i < fmpz_mat_nrows(A); i++)
    {
        for (j = 0; j < fmpz_mat_ncols(A); j++)
        {
            fmpz_sub_ui(fmpz_mat_entry(B, i, j),
                        fmpz_mat_entry(A, i, j), x);
        }
    }
}
