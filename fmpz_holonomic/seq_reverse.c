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
fmpz_holonomic_seq_reverse(fmpz_holonomic_t res, const fmpz_holonomic_t op)
{
    long i, j, r;

    r = fmpz_holonomic_order(op);

    fmpz_holonomic_set(res, op);

    if (r > 0)
    {
        fmpz c = r;

        for (i = 0; i <= r; i++)
            for (j = 1; j < fmpz_poly_length(res->coeffs + i); j += 2)
                fmpz_neg((res->coeffs + i)->coeffs + j,
                         (res->coeffs + i)->coeffs + j);

        for (i = 0; i <= r; i++)
            fmpz_poly_taylor_shift(res->coeffs + i, res->coeffs + i, &c);

        for (i = 0; i <= r / 2; i++)
            fmpz_poly_swap(res->coeffs + i, res->coeffs + r - i);

        fmpz_holonomic_normalise_sign(res);
    }
}

