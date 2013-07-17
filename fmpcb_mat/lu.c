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
fmpcb_mat_lu(long * P, fmpcb_mat_t LU, const fmpcb_mat_t A, long prec)
{
    fmpcb_t d, e;
    fmpcb_ptr * a;
    long i, j, m, n, r, row, col;
    int result;

    m = fmpcb_mat_nrows(A);
    n = fmpcb_mat_ncols(A);

    result = 1;

    if (m == 0 || n == 0)
        return result;

    fmpcb_mat_set(LU, A);

    a = LU->rows;

    row = col = 0;
    for (i = 0; i < m; i++)
        P[i] = i;

    fmpcb_init(d);
    fmpcb_init(e);

    while (row < m && col < n)
    {
        r = fmpcb_mat_find_pivot_partial(LU, row, m, col);

        if (r == -1)
        {
            result = 0;
            break;
        }
        else if (r != row)
            fmpcb_mat_swap_rows(LU, P, row, r);

        fmpcb_set(d, a[row] + col);

        for (j = row + 1; j < m; j++)
        {
            fmpcb_div(e, a[j] + col, d, prec);
            fmpcb_neg(e, e);
            _fmpcb_vec_scalar_addmul(a[j] + col,
                a[row] + col, n - col, e, prec);
            fmpcb_zero(a[j] + col);
            fmpcb_neg(a[j] + row, e);
        }

        row++;
        col++;
    }

    fmpcb_clear(d);
    fmpcb_clear(e);

    return result;
}
