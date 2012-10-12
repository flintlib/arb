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
fmpz_holonomic_fun_set_pow_fmpq(fmpz_holonomic_t op, const fmpq_t e)
{
    if (fmpz_is_one(fmpq_denref(e)))
    {
        fmpz_holonomic_fun_set_pow_fmpz(op, fmpq_numref(e));
    }
    else
    {
        fmpz_holonomic_fit_length(op, 2);

        fmpz_poly_set_fmpz(op->coeffs + 0, fmpq_numref(e));
        fmpz_poly_neg(op->coeffs + 0, op->coeffs + 0);

        fmpz_poly_zero(op->coeffs + 1);
        fmpz_poly_set_coeff_fmpz(op->coeffs + 1, 1, fmpq_denref(e));

        _fmpz_holonomic_set_length(op, 2);
    }
}

