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

#include "fmprb_mat.h"

int
fmprb_mat_solve(fmprb_mat_t X, const fmprb_mat_t A, const fmprb_mat_t B, long prec)
{
    int result;
    long n, m, *perm;
    fmprb_mat_t LU;

    n = fmprb_mat_nrows(A);
    m = fmprb_mat_ncols(X);

    if (n == 0 || m == 0)
        return 1;

    perm = _perm_init(n);
    fmprb_mat_init(LU, n, n);

    result = fmprb_mat_lu(perm, LU, A, prec);

    if (result)
        fmprb_mat_solve_lu_precomp(X, perm, LU, B, prec);

    fmprb_mat_clear(LU);
    _perm_clear(perm);

    return result;
}
