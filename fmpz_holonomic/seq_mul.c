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
fmpz_holonomic_seq_mul(fmpz_holonomic_t res, const fmpz_holonomic_t op1, const fmpz_holonomic_t op2)
{
    if (op1->length <= 1 || op2->length <= 1)
    {
        fmpz_holonomic_one(res);
    }
    else if (op1->length != 2 || op2->length != 2)
    {
        printf("not implemented: product of non-hypergeometric sequences\n");
        abort();
    }
    else
    {
        fmpz_holonomic_fit_length(res, 2);

        /* todo: better to compute gcd first? */
        fmpz_poly_mul(res->coeffs + 0, op1->coeffs + 0, op2->coeffs + 0);
        fmpz_poly_mul(res->coeffs + 1, op1->coeffs + 1, op2->coeffs + 1);
        fmpz_poly_neg(res->coeffs + 1, res->coeffs + 1);
        _fmpz_holonomic_set_length(res, 2);

        fmpz_holonomic_seq_normalise(res);
    }
}

