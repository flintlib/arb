/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "arb_poly.h"

void
arb_poly_derivative(arb_poly_t z, const arb_poly_t x)
{
    long i, len;

    len = x->length;

    if (len <= 1)
    {
        arb_poly_zero(z);
        return;
    }

    _arb_poly_fit_length(z, len - 1);

    for (i = 1; i < len; i++)
        fmpz_mul_ui(arb_poly_coeffs(z) + i - 1, arb_poly_coeffs(x) + i, i);

    z->length = len - 1;

    fmpz_mul_ui(arb_poly_radref(z), arb_poly_radref(x), len - 1);
    fmpz_set(arb_poly_expref(z), arb_poly_expref(x));
}
