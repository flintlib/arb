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
_arb_poly_normalise(fmpz * coeffs, fmpz_t exp, fmpz_t rad, long len, long bits)
{
    if (len != 0)
    {
        long b, shift;

        b = _fmpz_vec_max_bits(coeffs, len);
        b = FLINT_ABS(b);
        shift = FLINT_MAX(b - bits, 0);

        if (shift > 0)
        {
            _fmpz_vec_scalar_tdiv_q_2exp(coeffs, coeffs, len, shift);
            fmpz_cdiv_q_2exp(rad, rad, shift);  /* must round up */
            fmpz_add_ui(rad, rad, 1UL);         /* error from rounding */
            fmpz_add_ui(exp, exp, shift);
        }
    }
}
