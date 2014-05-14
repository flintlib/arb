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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "arb_mat.h"

void
arb_mat_set_fmpz_mat(arb_mat_t dest, const fmpz_mat_t src)
{
    long i, j;

    if (arb_mat_ncols(dest) != 0)
    {
        for (i = 0; i < arb_mat_nrows(dest); i++)
            for (j = 0; j < arb_mat_ncols(dest); j++)
                arb_set_fmpz(arb_mat_entry(dest, i, j),
                    fmpz_mat_entry(src, i, j));
    }
}
