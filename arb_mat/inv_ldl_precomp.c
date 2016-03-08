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

void
arb_mat_inv_ldl_precomp(arb_mat_t X, const arb_mat_t L, slong prec)
{
    slong n;

    if (arb_mat_nrows(X) != arb_mat_nrows(L) ||
        arb_mat_ncols(X) != arb_mat_ncols(L))
    {
        flint_printf("arb_mat_inv_ldl_precomp: incompatible dimensions\n");
        abort();
    }

    if (arb_mat_is_empty(L))
        return;

    n = arb_mat_nrows(L);

    if (n == 1)
    {
        arb_inv(arb_mat_entry(X, 0, 0), arb_mat_entry(L, 0, 0), prec);
        return;
    }

    if (X == L)
    {
        flint_printf("arb_mat_inv_ldl_precomp: unsupported aliasing\n");
        abort();
    }

    /* invert a 2x2 or larger matrix given its L * D * L^T decomposition */
    {
        slong i, j, k;
        arb_struct *s;

        s = _arb_vec_init(n);
        for (i = 0; i < n; i++)
        {
            arb_inv(s + i, arb_mat_entry(L, i, i), prec);
        }

        arb_mat_zero(X);
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
                arb_set(arb_mat_entry(X, j, i),
                        arb_mat_entry(X, i, j));
            }
        }
        _arb_vec_clear(s, n);
    }
}
