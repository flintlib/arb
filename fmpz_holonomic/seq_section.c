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

#include "fmpz_holonomic.h"
#include "fmpz_poly_mat.h"

void
fmpz_holonomic_seq_section(fmpz_holonomic_t res, const fmpz_holonomic_t op, long m)
{
    long i, j, r, rows, cols, left_cols, right_cols, pos;
    fmpz_poly_mat_t mat;
    fmpz_poly_t t;

    r = fmpz_holonomic_order(op);

    if (m == 0)
    {
        if (r < 1)
            fmpz_holonomic_one(res);
        else
            fmpz_holonomic_seq_set_const(res);
        return;
    }
    else if (m == 1)
    {
        fmpz_holonomic_set(res, op);
        return;
    }
    else if (m < 0)
    {
        fmpz_holonomic_seq_section(res, op, -m);
        fmpz_holonomic_seq_reverse(res, res);
        return;
    }

    /*
    Algorithm: by shifting the input recurrence, we obtain a set of
    linear relations over c[m*k], c[m*k+1], c[m*k+2]..., c[m*(k+r)].
    We use Gaussian elimination to eliminate the offsets that
    are not multiples of m, obtaining a relation for c[m*k], c[m*(k+1)],
    ..., c[m*(k+r)].
    */

    rows = m * r - r + 1;
    cols = m * r + 1;

    /* we put the non-multiples of m on the left of the matrix
       and the multiples on the right */
    right_cols = r + 1;
    left_cols = cols - right_cols;

    fmpz_poly_mat_init(mat, rows, cols);
    fmpz_poly_init(t);

    for (i = 0; i < rows; i++)
    {
        fmpz_poly_zero(t);
        fmpz_poly_set_coeff_si(t, 0, i);
        fmpz_poly_set_coeff_si(t, 1, m);

        for (j = 0; j < cols; j++)
        {
            /* zero */
            if (j < i || j - i > r)
                continue;

            /* insert in the right or the left part of the matrix */
            if (j % m == 0)
                pos = left_cols + (j / m);
            else
                pos = j - (j / m) - 1;

            fmpz_poly_compose(mat->rows[i] + pos, op->coeffs + j - i, t);
        }
    }

    /* Do Gaussian elimination */
    fmpz_poly_mat_fflu(mat, t, NULL, mat, 0);

    /* The relation is the last row */
    fmpz_holonomic_fit_length(res, r + 1);
    for (i = 0; i <= r; i++)
        fmpz_poly_neg(res->coeffs + i, mat->rows[rows - 1] + left_cols + i);
    _fmpz_holonomic_set_length(res, r + 1);
    fmpz_holonomic_seq_normalise(res);

    fmpz_poly_mat_clear(mat);
    fmpz_poly_clear(t);
}

