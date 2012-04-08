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
arb_poly_integral(arb_poly_t z, const arb_poly_t x)
{
    long i, len, bits, shift;

    len = x->length;

    if (len == 0)
    {
        arb_poly_zero(z);
        return;
    }

    _arb_poly_fit_length(z, len + 1);

    bits = _fmpz_vec_max_bits(arb_poly_coeffs(x), len);
    bits = FLINT_ABS(bits);
    shift = arb_poly_prec(z) - bits + FLINT_BIT_COUNT(len);

    if (shift >= 0)
    {
        for (i = len; i > 0; i--)
        {
            fmpz_mul_2exp(arb_poly_coeffs(z) + i, arb_poly_coeffs(x) + i - 1, shift);
            fmpz_tdiv_q_ui(arb_poly_coeffs(z) + i, arb_poly_coeffs(z) + i, i);
        }

        fmpz_mul_2exp(arb_poly_radref(z), arb_poly_radref(x), shift);
        fmpz_add_ui(arb_poly_radref(z), arb_poly_radref(z), 1);
        fmpz_sub_ui(arb_poly_expref(z), arb_poly_expref(x), shift);
    }
    else
    {
        /* todo: shift down to match target precision? */
        for (i = len; i > 0; i--)
            fmpz_tdiv_q_ui(arb_poly_coeffs(z) + i, arb_poly_coeffs(x) + i - 1, i);

        fmpz_add_ui(arb_poly_radref(z), arb_poly_radref(x), 1);
        fmpz_set(arb_poly_expref(z), arb_poly_expref(x));
    }

    fmpz_zero(arb_poly_coeffs(z));

    z->length = len + 1;
}
