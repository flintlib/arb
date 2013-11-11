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

#include "fmpcb_mat.h"

int
fmpcb_mat_contains(const fmpcb_mat_t mat1, const fmpcb_mat_t mat2)
{
    long i, j;

    if ((fmpcb_mat_nrows(mat1) != fmpcb_mat_nrows(mat2)) ||
        (fmpcb_mat_ncols(mat1) != fmpcb_mat_ncols(mat2)))
        return 0;

    for (i = 0; i < fmpcb_mat_nrows(mat1); i++)
        for (j = 0; j < fmpcb_mat_ncols(mat1); j++)
            if (!fmpcb_contains(fmpcb_mat_entry(mat1, i, j), fmpcb_mat_entry(mat2, i, j)))
                return 0;

    return 1;
}
