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
fmpcb_mat_inv(fmpcb_mat_t X, const fmpcb_mat_t A, long prec)
{
    if (X == A)
    {
        int r;
        fmpcb_mat_t T;
        fmpcb_mat_init(T, fmpcb_mat_nrows(A), fmpcb_mat_ncols(A));
        r = fmpcb_mat_inv(T, A, prec);
        fmpcb_mat_swap(T, X);
        fmpcb_mat_clear(T);
        return r;
    }

    fmpcb_mat_one(X);
    return fmpcb_mat_solve(X, A, X, prec);
}
