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

#include "acb_mat.h"

void
acb_mat_bound_frobenius_norm(mag_t b, const acb_mat_t A)
{
    slong i, j, r, c;
    mag_t t;

    r = acb_mat_nrows(A);
    c = acb_mat_ncols(A);

    mag_zero(b);

    if (r == 0 || c == 0)
        return;

    mag_init(t);

    for (i = 0; i < r; i++)
    {
        for (j = 0; j < c; j++)
        {
            acb_ptr z = acb_mat_entry(A, i, j);
            arb_get_mag(t, acb_realref(z));
            mag_addmul(b, t, t);
            arb_get_mag(t, acb_imagref(z));
            mag_addmul(b, t, t);
        }
    }

    mag_sqrt(b, b);

    mag_clear(t);
}
