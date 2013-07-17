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

long
fmpcb_mat_gauss_partial(fmpcb_mat_t A, long prec)
{
    fmpcb_t e;
    fmpcb_ptr * a;
    long j, m, n, r, rank, row, col, sign;

    m = A->r;
    n = A->c;
    a = A->rows;
    rank = row = col = 0;
    sign = 1;

    fmpcb_init(e);

    while (row < m && col < n)
    {
        r = fmpcb_mat_find_pivot_partial(A, row, m, col);

        if (r == -1)
        {
            break;
        }
        else if (r != row)
        {
            fmpcb_mat_swap_rows(A, NULL, row, r);
            sign *= -1;
        }

        rank++;

        for (j = row + 1; j < m; j++)
        {
            fmpcb_div(e, a[j] + col, a[row] + col, prec);
            fmpcb_neg(e, e);
            _fmpcb_vec_scalar_addmul(a[j] + col + 1, a[row] + col + 1, n - col - 1, e, prec);
        }

        row++;
        col++;
    }

    fmpcb_clear(e);

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
fmpcb_vec_get_fmpr_2norm_squared_bound(fmpr_t s, fmpcb_srcptr vec, long len, long prec)
{
    long i;
    fmpr_t t;

    fmpr_init(t);
    fmpr_zero(s);

    for (i = 0; i < len; i++)
    {
        fmprb_get_fmpr_abs_ubound(t, fmpcb_realref(vec + i), prec);
        fmpr_addmul(s, t, t, prec, FMPR_RND_UP);
        fmprb_get_fmpr_abs_ubound(t, fmpcb_imagref(vec + i), prec);
        fmpr_addmul(s, t, t, prec, FMPR_RND_UP);
    }

    fmpr_clear(t);
}

void
fmpcb_mat_det_inplace(fmpcb_t det, fmpcb_mat_t A, long prec)
{
    long i, n, sign, rank;

    n = fmpcb_mat_nrows(A);
    rank = fmpcb_mat_gauss_partial(A, prec);
    sign = (rank < 0) ? -1 : 1;
    rank = FLINT_ABS(rank);

    fmpcb_set_si(det, sign);
    for (i = 0; i < rank; i++)
        fmpcb_mul(det, det, fmpcb_mat_entry(A, i, i), prec);

    /* bound unreduced part using Hadamard's inequality */
    if (rank < n)
    {
        fmpr_t t;
        fmprb_t d;
        fmpcb_t e;

        fmpr_init(t);
        fmprb_init(d);
        fmpcb_init(e);

        fmpr_one(fmprb_radref(d));

        for (i = rank; i < n; i++)
        {
            fmpcb_vec_get_fmpr_2norm_squared_bound(t, A->rows[i] + rank, 
                n - rank, FMPRB_RAD_PREC);
            fmpr_mul(fmprb_radref(d), fmprb_radref(d), t, FMPRB_RAD_PREC, FMPR_RND_UP);
        }

        /* now d contains the absolute value of the determinant */
        fmpr_sqrt(fmprb_radref(d), fmprb_radref(d), FMPRB_RAD_PREC, FMPR_RND_UP);

        /* multiply by interval containing the unit disc */
        fmpr_set_ui(fmprb_radref(fmpcb_realref(e)), 1);
        fmpr_set_ui(fmprb_radref(fmpcb_imagref(e)), 1);
        fmpcb_mul_fmprb(e, e, d, prec);

        fmpcb_mul(det, det, e, prec);

        fmpcb_clear(e);
        fmprb_clear(d);
        fmpr_clear(t);
    }
}

void
fmpcb_mat_det(fmpcb_t det, const fmpcb_mat_t A, long prec)
{
    long n = fmpcb_mat_nrows(A);

    if (n == 0)
    {
        fmpcb_one(det);
    }
    else if (n == 1)
    {
        fmpcb_set(det, fmpcb_mat_entry(A, 0, 0));
    }
    else if (n == 2)
    {
        fmpcb_mul(det, fmpcb_mat_entry(A, 0, 0), fmpcb_mat_entry(A, 1, 1), prec);
        fmpcb_submul(det, fmpcb_mat_entry(A, 0, 1), fmpcb_mat_entry(A, 1, 0), prec);
    }
    else
    {
        fmpcb_mat_t T;
        fmpcb_mat_init(T, fmpcb_mat_nrows(A), fmpcb_mat_ncols(A));
        fmpcb_mat_set(T, A);
        fmpcb_mat_det_inplace(det, T, prec);
        fmpcb_mat_clear(T);
    }
}
