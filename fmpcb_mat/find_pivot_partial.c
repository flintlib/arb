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
fmpcb_cmpabs_approx(const fmpcb_t x, const fmpcb_t y)
{
    const fmpr_struct *xm, *ym;

    if (fmpr_cmpabs(fmprb_midref(fmpcb_realref(x)), fmprb_midref(fmpcb_imagref(x))) >= 0)
        xm = fmprb_midref(fmpcb_realref(x));
    else
        xm = fmprb_midref(fmpcb_imagref(x));

    if (fmpr_cmpabs(fmprb_midref(fmpcb_realref(y)), fmprb_midref(fmpcb_imagref(y))) >= 0)
        ym = fmprb_midref(fmpcb_realref(y));
    else
        ym = fmprb_midref(fmpcb_imagref(y));

    return fmpr_cmpabs(xm, ym);
}

long
fmpcb_mat_find_pivot_partial(const fmpcb_mat_t mat,
                                    long start_row, long end_row, long c)
{
    long best_row, i;

    best_row = -1;

    for (i = start_row; i < end_row; i++)
    {
        if (!fmpcb_contains_zero(fmpcb_mat_entry(mat, i, c)))
        {
            if (best_row == -1)
            {
                best_row = i;
            }
            /* todo: should take the radius into account */
            else if (fmpcb_cmpabs_approx(fmpcb_mat_entry(mat, i, c), fmpcb_mat_entry(mat, best_row, c)) > 0)
            {
                best_row = i;
            }
        }
    }
    return best_row;
}

