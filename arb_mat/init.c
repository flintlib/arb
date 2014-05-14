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
arb_mat_init(arb_mat_t mat, long r, long c)
{
    if (r != 0 && c != 0)
    {
        long i;

        mat->entries = _arb_vec_init(r * c);
        mat->rows = (arb_ptr *) flint_malloc(r * sizeof(arb_ptr));

        for (i = 0; i < r; i++)
            mat->rows[i] = mat->entries + i * c;
    }
    else
        mat->entries = NULL;

    mat->r = r;
    mat->c = c;
}
