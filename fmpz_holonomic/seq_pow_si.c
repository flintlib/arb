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
fmpz_holonomic_seq_pow_si(fmpz_holonomic_t res, const fmpz_holonomic_t op, long e)
{
    if (e == 0)
    {
        fmpz_holonomic_seq_set_const(res);
    }
    else if (op->length == 0 || op->length == 1 || e == 1)
    {
        fmpz_holonomic_set(res, op);
    }
    else if (op->length != 2)
    {
        printf("not implemented: power of non-hypergeometric sequence\n");
        abort();
    }
    else if (e < 0)
    {
        fmpz_holonomic_seq_pow_si(res, op, -e);
        fmpz_poly_swap(res->coeffs, res->coeffs + 1);
        fmpz_holonomic_seq_normalise(res);  /* signs only? */
    }
    else
    {
        fmpz_holonomic_fit_length(res, 2);
        fmpz_poly_pow(res->coeffs + 0, op->coeffs + 0, e);
        fmpz_poly_pow(res->coeffs + 1, op->coeffs + 1, e);

        if (e % 2 == 0)
            fmpz_poly_neg(res->coeffs + 0, res->coeffs + 0);
        _fmpz_holonomic_set_length(res, 2);
    }
}

