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
arb_mat_solve(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, long prec)
{
    int result;
    long n, m, *perm;
    arb_mat_t LU;

    n = arb_mat_nrows(A);
    m = arb_mat_ncols(X);

    if (n == 0 || m == 0)
        return 1;

    perm = _perm_init(n);
    arb_mat_init(LU, n, n);

    result = arb_mat_lu(perm, LU, A, prec);

    if (result)
        arb_mat_solve_lu_precomp(X, perm, LU, B, prec);

    arb_mat_clear(LU);
    _perm_clear(perm);

    return result;
}
