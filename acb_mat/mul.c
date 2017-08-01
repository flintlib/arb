/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

void
acb_mat_mul(acb_mat_t C, const acb_mat_t A, const acb_mat_t B, slong prec)
{
    if (flint_get_num_threads() > 1 &&
        ((double) acb_mat_nrows(A) *
         (double) acb_mat_nrows(B) *
         (double) acb_mat_ncols(B) *
         (double) prec > 100000))
    {
        acb_mat_mul_threaded(C, A, B, prec);
    }
    else
    {
        acb_mat_mul_classical(C, A, B, prec);
    }
}

