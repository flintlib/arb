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

#include "arb_mat.h"

void
arb_mat_mul(arb_mat_t C, const arb_mat_t A, const arb_mat_t B, long prec)
{
    if (flint_get_num_threads() > 1 &&
        ((double) arb_mat_nrows(A) *
         (double) arb_mat_nrows(B) *
         (double) arb_mat_ncols(B) *
         (double) prec > 100000))
    {
        arb_mat_mul_threaded(C, A, B, prec);
    }
    else
    {
        arb_mat_mul_classical(C, A, B, prec);
    }
}

