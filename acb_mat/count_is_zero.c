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

#include "acb_mat.h"

slong
acb_mat_count_is_zero(const acb_mat_t mat)
{
    slong nz;
    slong i, j;
    nz = 0;
    for (i = 0; i < acb_mat_nrows(mat); i++)
    {
        for (j = 0; j < acb_mat_ncols(mat); j++)
        {
            if (acb_is_zero(acb_mat_entry(mat, i, j)))
            {
                nz++;
            }
        }
    }
    return nz;
}
