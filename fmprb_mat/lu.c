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

static __inline__ void
fmprb_mat_swap_rows(fmprb_mat_t mat, long * perm, long r, long s)
{
    if (r != s)
    {
        fmprb_struct * u;
        long t;

        if (perm != NULL)
        {
            t = perm[s];
            perm[s] = perm[r];
            perm[r] = t;
        }

        u = mat->rows[s];
        mat->rows[s] = mat->rows[r];
        mat->rows[r] = u;
    }
}

static __inline__ long
fmprb_mat_find_pivot_partial(const fmprb_mat_t mat,
                                    long start_row, long end_row, long c)
{
    long best_row, i;

    best_row = -1;

    for (i = start_row; i < end_row; i++)
    {
        if (!fmprb_contains_zero(fmprb_mat_entry(mat, i, c)))
        {
            if (best_row == -1)
            {
                best_row = i;
            }
            /* todo: should take the radius into account */
            else if (fmpr_cmpabs(fmprb_midref(fmprb_mat_entry(mat, i, c)),
                    fmprb_midref(fmprb_mat_entry(mat, best_row, c))) > 0)
            {
                best_row = i;
            }
        }
    }

    return best_row;
}

int
fmprb_mat_lu(long * P, fmprb_mat_t LU, const fmprb_mat_t A, long prec)
{
    fmprb_t d, e;
    fmprb_struct ** a;
    long i, j, m, n, r, length, row, col;
    int result;

    m = fmprb_mat_nrows(A);
    n = fmprb_mat_ncols(A);

    result = 1;

    if (m == 0 || n == 0)
        return result;

    fmprb_mat_set(LU, A);

    a = LU->rows;

    row = col = 0;
    for (i = 0; i < m; i++)
        P[i] = i;

    fmprb_init(d);
    fmprb_init(e);

    while (row < m && col < n)
    {
        r = fmprb_mat_find_pivot_partial(LU, row, m, col);

        if (r == -1)
        {
            result = 0;
            break;
        }
        else if (r != row)
            fmprb_mat_swap_rows(LU, P, row, r);

        fmprb_set(d, a[row] + col);

        for (j = row + 1; j < m; j++)
        {
            fmprb_div(e, a[j] + col, d, prec);
            fmprb_neg(e, e);
            _fmprb_vec_scalar_addmul(a[j] + col,
                a[row] + col, n - col, e, prec);
            fmprb_zero(a[j] + col);
            fmprb_neg(a[j] + row, e);
        }

        row++;
        col++;
    }

    fmprb_clear(d);
    fmprb_clear(e);

    return result;
}
