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

    Copyright (C) 2015 Arb authors

******************************************************************************/

#include "acb_mat.h"

void
acb_mat_trace(acb_t trace, const acb_mat_t mat, slong prec)
{
    slong i, n = acb_mat_nrows(mat);

    if (n == 0)
    {
        acb_zero(trace);
    }
    else
    {
        acb_set(trace, acb_mat_entry(mat, 0, 0));
        for (i = 1; i < n; i++)
        {
            acb_add(trace, trace, acb_mat_entry(mat, i, i), prec);
        }
    }
}
