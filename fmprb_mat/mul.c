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

void
fmprb_mat_mul(fmprb_mat_t C, const fmprb_mat_t A, const fmprb_mat_t B, long prec)
{
    long ar, ac, br, bc, i, j, k;

    ar = fmprb_mat_nrows(A);
    ac = fmprb_mat_ncols(A);
    br = fmprb_mat_nrows(B);
    bc = fmprb_mat_ncols(B);

    if (ac != br || ar != fmprb_mat_nrows(C) || bc != fmprb_mat_ncols(C))
    {
        printf("fmprb_mat_mul: incompatible dimensions\n");
        abort();
    }

    if (br == 0)
    {
        fmprb_mat_zero(C);
        return;
    }

    if (A == C || B == C)
    {
        fmprb_mat_t T;
        fmprb_mat_init(T, ar, bc);
        fmprb_mat_mul(T, A, B, prec);
        fmprb_mat_set(C, T);
        fmprb_mat_clear(T);
        return;
    }

    for (i = 0; i < ar; i++)
    {
        for (j = 0; j < bc; j++)
        {
            fmprb_mul(fmprb_mat_entry(C, i, j),
                      fmprb_mat_entry(A, i, 0),
                      fmprb_mat_entry(B, 0, j), prec);

            for (k = 1; k < br; k++)
            {
                fmprb_addmul(fmprb_mat_entry(C, i, j),
                             fmprb_mat_entry(A, i, k),
                             fmprb_mat_entry(B, k, j), prec);
            }
        }
    }
}
