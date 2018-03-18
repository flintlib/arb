/*
    Copyright (C) 2018 arbguest

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

void
_arb_vec_approx_scalar_addmul(arb_ptr res, arb_srcptr vec,
    slong len, const arb_t c, slong prec)
{
    slong i;
    for (i = 0; i < len; i++)
        arf_addmul(arb_midref(res + i),
                   arb_midref(vec + i), arb_midref(c), prec, ARF_RND_DOWN);
}

int
arb_mat_approx_lu(slong * P, arb_mat_t LU, const arb_mat_t A, slong prec)
{
    arf_t d;
    arb_t e;
    arb_ptr * a;
    slong i, j, m, n, r, row, col;
    int result;

    if (arb_mat_is_empty(A))
        return 1;

    m = arb_mat_nrows(A);
    n = arb_mat_ncols(A);

    arb_mat_get_mid(LU, A);

    a = LU->rows;

    row = col = 0;
    for (i = 0; i < m; i++)
        P[i] = i;

    arf_init(d);
    arb_init(e);

    result = 1;

    while (row < m && col < n)
    {
        r = arb_mat_find_pivot_partial(LU, row, m, col);

        if (r == -1)
        {
            result = 0;
            break;
        }
        else if (r != row)
            arb_mat_swap_rows(LU, P, row, r);

        arf_set(d, arb_midref(a[row] + col));

        for (j = row + 1; j < m; j++)
        {
            arf_div(arb_midref(e), arb_midref(a[j] + col), d, prec, ARB_RND);
            arb_neg(e, e);
            _arb_vec_approx_scalar_addmul(a[j] + col,
                a[row] + col, n - col, e, prec);
            arf_zero(arb_midref(a[j] + col));
            arb_neg(a[j] + row, e);
        }

        row++;
        col++;
    }

    arf_clear(d);
    arb_clear(e);

    return result;
}
