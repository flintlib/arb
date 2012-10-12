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
fmpz_holonomic_fun_set_asin_acos(fmpz_holonomic_t op)
{
    fmpz_holonomic_fit_length(op, 3);
    fmpz_poly_set_si(op->coeffs + 0, 0);
    fmpz_poly_set_si2(op->coeffs + 1, 0, 1);
    fmpz_poly_set_si3(op->coeffs + 2, -1, 0, 1);
    _fmpz_holonomic_set_length(op, 3);
}

