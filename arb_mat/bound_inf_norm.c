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

#include "arb_mat.h"

void
arb_mat_bound_inf_norm(arf_t b, const arb_mat_t A, long prec)
{
    long i, j, r, c;

    arf_t s, t;

    r = arb_mat_nrows(A);
    c = arb_mat_ncols(A);

    arf_zero(b);

    if (r == 0 || c == 0)
        return;

    arf_init(s);
    arf_init(t);

    for (i = 0; i < r; i++)
    {
        arf_zero(s);

        for (j = 0; j < c; j++)
        {
            arb_get_abs_ubound_arf(t, arb_mat_entry(A, i, j), prec);
            arf_add(s, s, t, prec, ARF_RND_UP);
        }

        arf_max(b, b, s);
    }

    arf_clear(s);
    arf_clear(t);
}

