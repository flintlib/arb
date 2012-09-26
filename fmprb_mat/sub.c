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

#include "fmprb_mat.h"

void
fmprb_mat_sub(fmprb_mat_t res,
        const fmprb_mat_t mat1, const fmprb_mat_t mat2, long prec)
{
    long i, j;

    for (i = 0; i < fmprb_mat_nrows(mat1); i++)
        for (j = 0; j < fmprb_mat_ncols(mat1); j++)
            fmprb_sub(fmprb_mat_entry(res, i, j),
                fmprb_mat_entry(mat1, i, j),
                fmprb_mat_entry(mat2, i, j), prec);
}
