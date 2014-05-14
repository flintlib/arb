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

long
arb_mat_find_pivot_partial(const arb_mat_t mat,
                                    long start_row, long end_row, long c)
{
    long best_row, i;

    best_row = -1;

    for (i = start_row; i < end_row; i++)
    {
        if (!arb_contains_zero(arb_mat_entry(mat, i, c)))
        {
            if (best_row == -1)
            {
                best_row = i;
            }
            /* todo: should take the radius into account */
            else if (arf_cmpabs(arb_midref(arb_mat_entry(mat, i, c)),
                                arb_midref(arb_mat_entry(mat, best_row, c))) > 0)
            {
                best_row = i;
            }
        }
    }

    return best_row;
}
