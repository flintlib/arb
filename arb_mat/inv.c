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

int
arb_mat_inv(arb_mat_t X, const arb_mat_t A, long prec)
{
    if (X == A)
    {
        int r;
        arb_mat_t T;
        arb_mat_init(T, arb_mat_nrows(A), arb_mat_ncols(A));
        r = arb_mat_inv(T, A, prec);
        arb_mat_swap(T, X);
        arb_mat_clear(T);
        return r;
    }

    arb_mat_one(X);
    return arb_mat_solve(X, A, X, prec);
}
