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

void
fmpz_holonomic_fit_length(fmpz_holonomic_t op, long len)
{
    long i;

    if (len > op->alloc)
    {
        if (len < 2 * op->alloc)
            len = 2 * op->alloc;

        op->coeffs = flint_realloc(op->coeffs,
            len * sizeof(fmpz_poly_struct));

        for (i = op->alloc; i < len; i++)
            fmpz_poly_init(op->coeffs + i);

        op->alloc = len;
    }
}

