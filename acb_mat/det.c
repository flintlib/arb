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

long
acb_mat_gauss_partial(acb_mat_t A, long prec)
{
    acb_t e;
    acb_ptr * a;
    long j, m, n, r, rank, row, col, sign;

    m = A->r;
    n = A->c;
    a = A->rows;
    rank = row = col = 0;
    sign = 1;

    acb_init(e);

    while (row < m && col < n)
    {
        r = acb_mat_find_pivot_partial(A, row, m, col);

        if (r == -1)
        {
            break;
        }
        else if (r != row)
        {
            acb_mat_swap_rows(A, NULL, row, r);
            sign *= -1;
        }

        rank++;

        for (j = row + 1; j < m; j++)
        {
            acb_div(e, a[j] + col, a[row] + col, prec);
            acb_neg(e, e);
            _acb_vec_scalar_addmul(a[j] + col + 1, a[row] + col + 1, n - col - 1, e, prec);
        }

        row++;
        col++;
    }

    acb_clear(e);

    return rank * sign;
}

void
acb_vec_get_arf_2norm_squared_bound(arf_t s, acb_srcptr vec, long len, long prec)
{
    long i;
    arf_t t;

    arf_init(t);
    arf_zero(s);

    for (i = 0; i < len; i++)
    {
        arb_get_abs_ubound_arf(t, acb_realref(vec + i), prec);
        arf_addmul(s, t, t, prec, ARF_RND_UP);
        arb_get_abs_ubound_arf(t, acb_imagref(vec + i), prec);
        arf_addmul(s, t, t, prec, ARF_RND_UP);
    }

    arf_clear(t);
}

void
acb_mat_det_inplace(acb_t det, acb_mat_t A, long prec)
{
    long i, n, sign, rank;
    int is_real;

    n = acb_mat_nrows(A);
    rank = acb_mat_gauss_partial(A, prec);
    sign = (rank < 0) ? -1 : 1;
    rank = FLINT_ABS(rank);

    acb_set_si(det, sign);
    for (i = 0; i < rank; i++)
        acb_mul(det, det, acb_mat_entry(A, i, i), prec);

    /* bound unreduced part using Hadamard's inequality */
    if (rank < n)
    {
        arf_t t;
        arf_t d;
        acb_t e;

        arf_init(t);
        arf_init(d);
        acb_init(e);

        arf_one(d);

        is_real = acb_mat_is_real(A);

        for (i = rank; i < n; i++)
        {
            acb_vec_get_arf_2norm_squared_bound(t, A->rows[i] + rank,  n - rank, MAG_BITS);
            arf_mul(d, d, t, MAG_BITS, ARF_RND_UP);
        }

        /* now d contains the absolute value of the determinant */
        arf_sqrt(d, d, MAG_BITS, ARF_RND_UP);

        /* multiply by disc with radius d */
        if (is_real)
        {
            arb_add_error_arf(acb_realref(e), d);
        }
        else
        {
            arb_add_error_arf(acb_realref(e), d);
            arb_add_error_arf(acb_imagref(e), d);
        }

        acb_mul(det, det, e, prec);

        acb_clear(e);
        arf_clear(d);
        arf_clear(t);
    }
}

void
acb_mat_det(acb_t det, const acb_mat_t A, long prec)
{
    long n = acb_mat_nrows(A);

    if (n == 0)
    {
        acb_one(det);
    }
    else if (n == 1)
    {
        acb_set(det, acb_mat_entry(A, 0, 0));
    }
    else if (n == 2)
    {
        acb_mul(det, acb_mat_entry(A, 0, 0), acb_mat_entry(A, 1, 1), prec);
        acb_submul(det, acb_mat_entry(A, 0, 1), acb_mat_entry(A, 1, 0), prec);
    }
    else
    {
        acb_mat_t T;
        acb_mat_init(T, acb_mat_nrows(A), acb_mat_ncols(A));
        acb_mat_set(T, A);
        acb_mat_det_inplace(det, T, prec);
        acb_mat_clear(T);
    }
}
