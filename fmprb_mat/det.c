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

long
fmprb_mat_gauss_partial(fmprb_mat_t A, long prec)
{
    fmprb_t e;
    fmprb_ptr * a;
    long j, m, n, r, rank, row, col, sign;

    m = A->r;
    n = A->c;
    a = A->rows;
    rank = row = col = 0;
    sign = 1;

    fmprb_init(e);

    while (row < m && col < n)
    {
        r = fmprb_mat_find_pivot_partial(A, row, m, col);

        if (r == -1)
        {
            break;
        }
        else if (r != row)
        {
            fmprb_mat_swap_rows(A, NULL, row, r);
            sign *= -1;
        }

        rank++;

        for (j = row + 1; j < m; j++)
        {
            fmprb_div(e, a[j] + col, a[row] + col, prec);
            fmprb_neg(e, e);
            _fmprb_vec_scalar_addmul(a[j] + col + 1, a[row] + col + 1, n - col - 1, e, prec);
        }

        row++;
        col++;
    }

    fmprb_clear(e);

    return rank * sign;
}

static __inline__ void
fmprb_get_fmpr_abs_ubound(fmpr_t u, const fmprb_t x, long prec)
{
    if (fmpr_sgn(fmprb_midref(x)) >= 0)
    {
        fmpr_add(u, fmprb_midref(x), fmprb_radref(x), prec, FMPR_RND_UP);
    }
    else
    {
        fmpr_sub(u, fmprb_midref(x), fmprb_radref(x), prec, FMPR_RND_UP);
        fmpr_neg(u, u);
    }
}

void
fmprb_vec_get_fmpr_2norm_squared_bound(fmpr_t s, fmprb_srcptr vec, long len, long prec)
{
    long i;
    fmpr_t t;

    fmpr_init(t);
    fmpr_zero(s);

    for (i = 0; i < len; i++)
    {
        fmprb_get_fmpr_abs_ubound(t, vec + i, prec);
        fmpr_addmul(s, t, t, prec, FMPR_RND_UP);
    }

    fmpr_clear(t);
}

void
fmprb_mat_det_inplace(fmprb_t det, fmprb_mat_t A, long prec)
{
    long i, n, sign, rank;

    n = fmprb_mat_nrows(A);
    rank = fmprb_mat_gauss_partial(A, prec);
    sign = (rank < 0) ? -1 : 1;
    rank = FLINT_ABS(rank);

    fmprb_set_si(det, sign);
    for (i = 0; i < rank; i++)
        fmprb_mul(det, det, fmprb_mat_entry(A, i, i), prec);

    /* bound unreduced part using Hadamard's inequality */
    if (rank < n)
    {
        fmpr_t t;
        fmprb_t d;

        fmpr_init(t);
        fmprb_init(d);

        fmpr_one(fmprb_radref(d));

        for (i = rank; i < n; i++)
        {
            fmprb_vec_get_fmpr_2norm_squared_bound(t, A->rows[i] + rank, 
                n - rank, FMPRB_RAD_PREC);
            fmpr_mul(fmprb_radref(d), fmprb_radref(d), t, FMPRB_RAD_PREC, FMPR_RND_UP);
        }

        fmpr_sqrt(fmprb_radref(d), fmprb_radref(d), FMPRB_RAD_PREC, FMPR_RND_UP);
        fmprb_mul(det, det, d, prec);

        fmprb_clear(d);
        fmpr_clear(t);
    }
}

void
fmprb_mat_det(fmprb_t det, const fmprb_mat_t A, long prec)
{
    long n = fmprb_mat_nrows(A);

    if (n == 0)
    {
        fmprb_one(det);
    }
    else if (n == 1)
    {
        fmprb_set(det, fmprb_mat_entry(A, 0, 0));
    }
    else if (n == 2)
    {
        fmprb_mul(det, fmprb_mat_entry(A, 0, 0), fmprb_mat_entry(A, 1, 1), prec);
        fmprb_submul(det, fmprb_mat_entry(A, 0, 1), fmprb_mat_entry(A, 1, 0), prec);
    }
    else
    {
        fmprb_mat_t T;
        fmprb_mat_init(T, fmprb_mat_nrows(A), fmprb_mat_ncols(A));
        fmprb_mat_set(T, A);
        fmprb_mat_det_inplace(det, T, prec);
        fmprb_mat_clear(T);
    }
}
