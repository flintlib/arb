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

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include "fmprb_mat.h"

void
fmprb_mat_bound_inf_norm(fmpr_t b, const fmprb_mat_t A, long prec)
{
    long i, j, r, c;

    fmpr_t s, t;

    r = fmprb_mat_nrows(A);
    c = fmprb_mat_ncols(A);

    fmpr_zero(b);

    if (r == 0 || c == 0)
        return;

    fmpr_init(s);
    fmpr_init(t);

    for (i = 0; i < r; i++)
    {
        fmpr_zero(s);

        for (j = 0; j < c; j++)
        {
            fmprb_get_abs_ubound_fmpr(t, fmprb_mat_entry(A, i, j), prec);
            fmpr_add(s, s, t, prec, FMPR_RND_UP);
        }

        fmpr_max(b, b, s);
    }

    fmpr_clear(s);
    fmpr_clear(t);
}

