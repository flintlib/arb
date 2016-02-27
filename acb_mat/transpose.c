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

    Copyright (C) 2015 Tommy Hofmann

******************************************************************************/

#include "acb_mat.h"

void
acb_mat_transpose(acb_mat_t B, const acb_mat_t A)
{
    slong i, j;

    if (acb_mat_nrows(B) != acb_mat_ncols(A) || acb_mat_ncols(B) != acb_mat_nrows(A))
    {
        flint_printf("Exception (acb_mat_transpose). Incompatible dimensions.\n");
        abort();
    }

    if (acb_mat_is_empty(A))
        return;

    if (A == B)  /* In-place, guaranteed to be square */
    {
        for (i = 0; i < acb_mat_nrows(A) - 1; i++)
        {
            for (j = i + 1; j < acb_mat_ncols(A); j++)
            {
                acb_swap(acb_mat_entry(A, i, j), acb_mat_entry(A, j, i));
            }
        }
    }
    else  /* Not aliased; general case */
    {
        for (i = 0; i < acb_mat_nrows(B); i++)
            for (j = 0; j < acb_mat_ncols(B); j++)
                acb_set(acb_mat_entry(B, i, j), acb_mat_entry(A, j, i));
    }
}

