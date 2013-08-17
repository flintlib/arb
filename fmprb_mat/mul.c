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
    if (flint_get_num_threads() > 1 &&
        ((double) fmprb_mat_nrows(A) *
         (double) fmprb_mat_nrows(B) *
         (double) fmprb_mat_ncols(B) *
         (double) prec > 100000))
    {
        fmprb_mat_mul_threaded(C, A, B, prec);
    }
    else
    {
        fmprb_mat_mul_classical(C, A, B, prec);
    }
}

