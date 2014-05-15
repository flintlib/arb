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
acb_mat_inv(acb_mat_t X, const acb_mat_t A, long prec)
{
    if (X == A)
    {
        int r;
        acb_mat_t T;
        acb_mat_init(T, acb_mat_nrows(A), acb_mat_ncols(A));
        r = acb_mat_inv(T, A, prec);
        acb_mat_swap(T, X);
        acb_mat_clear(T);
        return r;
    }

    acb_mat_one(X);
    return acb_mat_solve(X, A, X, prec);
}
