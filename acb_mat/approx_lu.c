/*
    Copyright (C) 2018 arbguest

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

static void
_acb_approx_mul(acb_t res, const acb_t x, const acb_t y, slong prec)
{
    arf_complex_mul(arb_midref(acb_realref(res)), arb_midref(acb_imagref(res)),
        arb_midref(acb_realref(x)), arb_midref(acb_imagref(x)), 
        arb_midref(acb_realref(y)), arb_midref(acb_imagref(y)), prec, ARB_RND);
}

static void
_acb_approx_inv(acb_t z, const acb_t x, slong prec)
{
    arf_set(arb_midref(acb_realref(z)), arb_midref(acb_realref(x)));
    arf_set(arb_midref(acb_imagref(z)), arb_midref(acb_imagref(x)));

    mag_zero(arb_radref(acb_realref(z)));
    mag_zero(arb_radref(acb_imagref(z)));

    acb_inv(z, z, prec);

    mag_zero(arb_radref(acb_realref(z)));
    mag_zero(arb_radref(acb_imagref(z)));
}

static void
_acb_vec_approx_scalar_addmul(acb_ptr res, acb_srcptr vec,
    slong len, const acb_t c, slong prec)
{
    acb_t t;
    slong i;
    acb_init(t);

    for (i = 0; i < len; i++)
    {
        _acb_approx_mul(t, vec + i, c, prec);

        arf_add(arb_midref(acb_realref(res + i)),
            arb_midref(acb_realref(res + i)), 
            arb_midref(acb_realref(t)), prec, ARB_RND);
        arf_add(arb_midref(acb_imagref(res + i)),
            arb_midref(acb_imagref(res + i)), 
            arb_midref(acb_imagref(t)), prec, ARB_RND);
    }

    acb_clear(t);
}

int
acb_mat_approx_lu(slong * P, acb_mat_t LU, const acb_mat_t A, slong prec)
{
    acb_t d, e;
    acb_ptr * a;
    slong i, j, m, n, r, row, col;
    int result;

    if (acb_mat_is_empty(A))
        return 1;

    m = acb_mat_nrows(A);
    n = acb_mat_ncols(A);

    acb_mat_get_mid(LU, A);

    a = LU->rows;

    row = col = 0;
    for (i = 0; i < m; i++)
        P[i] = i;

    acb_init(d);
    acb_init(e);

    result = 1;

    while (row < m && col < n)
    {
        r = acb_mat_find_pivot_partial(LU, row, m, col);

        if (r == -1)
        {
            result = 0;
            break;
        }
        else if (r != row)
            acb_mat_swap_rows(LU, P, row, r);

        _acb_approx_inv(d, a[row] + col, prec);

        for (j = row + 1; j < m; j++)
        {
            _acb_approx_mul(e, a[j] + col, d, prec);
            acb_neg(e, e);
            _acb_vec_approx_scalar_addmul(a[j] + col,
                a[row] + col, n - col, e, prec);
            acb_zero(a[j] + col);
            acb_neg(a[j] + row, e);
        }

        row++;
        col++;
    }

    acb_clear(d);
    acb_clear(e);

    return result;
}
