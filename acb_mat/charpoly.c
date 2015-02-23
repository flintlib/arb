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

    Copyright (C) 2012 Sebastian Pancratz

******************************************************************************/

#include "acb_mat.h"

void _acb_mat_charpoly(acb_ptr cp, const acb_mat_t mat, long prec)
{
    const long n = mat->r;

    if (n == 0)
    {
        acb_one(cp);
    }
    else if (n == 1)
    {
        acb_neg(cp + 0, acb_mat_entry(mat, 0, 0));
        acb_one(cp + 1);
    }
    else
    {
        long i, j, k, t;
        acb_ptr a, A, s;

        a = _acb_vec_init(n * n);
        A = a + (n - 1) * n;

        _acb_vec_zero(cp, n + 1);
        acb_neg(cp + 0, acb_mat_entry(mat, 0, 0));

        for (t = 1; t < n; t++)
        {
            for (i = 0; i <= t; i++)
            {
                acb_set(a + 0 * n + i, acb_mat_entry(mat, i, t));
            }

            acb_set(A + 0, acb_mat_entry(mat, t, t));

            for (k = 1; k < t; k++)
            {
                for (i = 0; i <= t; i++)
                {
                    s = a + k * n + i;
                    acb_zero(s);
                    for (j = 0; j <= t; j++)
                        acb_addmul(s, acb_mat_entry(mat, i, j), a + (k - 1) * n + j, prec);
                }

                acb_set(A + k, a + k * n + t);
            }

            acb_zero(A + t);
            for (j = 0; j <= t; j++)
                acb_addmul(A + t, acb_mat_entry(mat, t, j), a + (t - 1) * n + j, prec);

            for (k = 0; k <= t; k++)
            {
                for (j = 0; j < k; j++)
                    acb_submul(cp + k, A + j, cp + (k - j - 1), prec);

                acb_sub(cp + k, cp + k, A + k, prec);
            }
        }

        /* Shift all coefficients up by one */
        for (i = n; i > 0; i--)
            acb_swap(cp + i, cp + (i - 1));

        acb_one(cp + 0);
        _acb_poly_reverse(cp, cp, n + 1, n + 1);
        _acb_vec_clear(a, n * n);
    }
}

void acb_mat_charpoly(acb_poly_t cp, const acb_mat_t mat, long prec)
{
    if (mat->r != mat->c)
    {
        flint_printf("Exception (acb_mat_charpoly).  Non-square matrix.\n");
        abort();
    }

    acb_poly_fit_length(cp, mat->r + 1);
    _acb_poly_set_length(cp, mat->r + 1);
    _acb_mat_charpoly(cp->coeffs, mat, prec);
}

