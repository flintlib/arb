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
fmpz_holonomic_normalise_sign(fmpz_holonomic_t op)
{
    if (op->length > 0)
    {
        long i;
        fmpz_poly_struct * poly = op->coeffs + op->length - 1;

        if (fmpz_sgn(poly->coeffs + poly->length - 1) < 0)
        {
            for (i = 0; i < op->length; i++)
                fmpz_poly_neg(op->coeffs + i, op->coeffs + i);
        }
    }
}

