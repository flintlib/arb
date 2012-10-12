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
fmpz_holonomic_seq_normalise_trailing(fmpz_holonomic_t op)
{
    long i, j;

    i = 0;
    while (i < op->length && fmpz_poly_is_zero(op->coeffs + i))
        i++;

    if (i > 0)
    {
        fmpz_t t;

        for (j = 0; j < op->length - i; j++)
            fmpz_poly_swap(op->coeffs + j, op->coeffs + j + i);

        op->length -= i;

        fmpz_init(t);
        fmpz_set_si(t, -i);
        for (j = 0; j < op->length; j++)
            fmpz_poly_taylor_shift(op->coeffs + j, op->coeffs + j, t);
        fmpz_clear(t);
    }
}

