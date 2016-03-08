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

    Copyright (C) 2016 Arb authors

******************************************************************************/

#include "arb_mat.h"

int
arb_mat_spd_inv(arb_mat_t X, const arb_mat_t A, slong prec)
{
    slong n;
    arb_mat_t L;
    int result;

    if (!arb_mat_is_square(A))
    {
        flint_printf("arb_mat_spd_inv: a square matrix is required\n");
        abort();
    }

    if (arb_mat_nrows(X) != arb_mat_nrows(A) ||
        arb_mat_ncols(X) != arb_mat_ncols(A))
    {
        flint_printf("arb_mat_spd_inv: incompatible dimensions\n");
        abort();
    }

    if (arb_mat_is_empty(A))
        return 1;

    n = arb_mat_nrows(A);

    if (n == 1)
    {
        if (arb_is_positive(arb_mat_entry(A, 0, 0)))
        {
            arb_inv(arb_mat_entry(X, 0, 0), arb_mat_entry(A, 0, 0), prec);
            return 1;
        }
        else
        {
            return 0;
        }
    }

    arb_mat_init(L, n, n);
    arb_mat_set(L, A);

    if (!_arb_mat_cholesky_banachiewicz(L, prec))
    {
        result = 0;
    }
    else
    {
        slong i, j, k;
        arb_struct *s;
        arb_mat_zero(X);
        s = _arb_vec_init(n);
        for (i = 0; i < n; i++)
        {
            arb_inv(s + i, arb_mat_entry(L, i, i), prec);
        }
        for (j = n-1; j >= 0; j--)
        {
            for (i = j; i >= 0; i--)
            {
                if (i == j)
                {
                    arb_set(arb_mat_entry(X, i, j), s + i);
                }
                else
                {
                    arb_zero(arb_mat_entry(X, i, j));
                }
                for (k = i + 1; k < n; k++)
                {
                    arb_submul(arb_mat_entry(X, i, j),
                               arb_mat_entry(L, k, i),
                               arb_mat_entry(X, k, j), prec);
                }
                arb_div(arb_mat_entry(X, i, j),
                        arb_mat_entry(X, i, j),
                        arb_mat_entry(L, i, i), prec);
                arb_set(arb_mat_entry(X, j, i),
                        arb_mat_entry(X, i, j));
            }
        }

        _arb_vec_clear(s, n);
        result = 1;
    }

    arb_mat_clear(L);
    return result;
}
