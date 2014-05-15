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

#include "acb_mat.h"

int
acb_cmpabs_approx(const acb_t x, const acb_t y)
{
    const arf_struct *xm, *ym;

    if (arf_cmpabs(arb_midref(acb_realref(x)), arb_midref(acb_imagref(x))) >= 0)
        xm = arb_midref(acb_realref(x));
    else
        xm = arb_midref(acb_imagref(x));

    if (arf_cmpabs(arb_midref(acb_realref(y)), arb_midref(acb_imagref(y))) >= 0)
        ym = arb_midref(acb_realref(y));
    else
        ym = arb_midref(acb_imagref(y));

    return arf_cmpabs(xm, ym);
}

long
acb_mat_find_pivot_partial(const acb_mat_t mat,
                                    long start_row, long end_row, long c)
{
    long best_row, i;

    best_row = -1;

    for (i = start_row; i < end_row; i++)
    {
        if (!acb_contains_zero(acb_mat_entry(mat, i, c)))
        {
            if (best_row == -1)
            {
                best_row = i;
            }
            /* todo: should take the radius into account */
            else if (acb_cmpabs_approx(acb_mat_entry(mat, i, c), acb_mat_entry(mat, best_row, c)) > 0)
            {
                best_row = i;
            }
        }
    }
    return best_row;
}

