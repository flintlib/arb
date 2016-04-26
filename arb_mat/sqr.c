/*
    Copyright (C) 2015 Arb authors

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

void
arb_mat_sqr(arb_mat_t B, const arb_mat_t A, slong prec)
{
    double d = (double) arb_mat_nrows(A);

    if (flint_get_num_threads() > 1 && d*d*d * (double) prec > 100000)
    {
        arb_mat_mul_threaded(B, A, A, prec);   /* todo: sqr threaded */
    }
    else
    {
        arb_mat_sqr_classical(B, A, prec);
    }
}

