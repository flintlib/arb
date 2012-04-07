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
arb_poly_set(arb_poly_t z, const arb_poly_t x)
{
    if (z == x)
        return;

    if (x->length == 0)
    {
        arb_poly_zero(z);
        return;
    }

    _arb_poly_fit_length(z, x->length);
    z->length = x->length;

    if (z->prec >= x->prec)
    {
        _fmpz_vec_set(z->coeffs, x->coeffs, x->length);
        fmpz_set(&z->exp, &x->exp);
        fmpz_set(&z->rad, &x->rad);
    }
    else
    {
        long bits;

        bits = _fmpz_vec_max_bits(x->coeffs, x->length);
        bits = FLINT_ABS(bits);

        if (bits <= z->prec)
        {
            _fmpz_vec_set(z->coeffs, x->coeffs, x->length);
            fmpz_set(&z->exp, &x->exp);
            fmpz_set(&z->rad, &x->rad);
        }
        else
        {
            _fmpz_vec_scalar_tdiv_q_2exp(z->coeffs, x->coeffs, x->length,
                bits - z->prec);
            fmpz_cdiv_q_2exp(&z->rad, &x->rad, bits - z->prec);
            fmpz_add_ui(&z->rad, &z->rad, 1);   /* rounding */
            fmpz_add_ui(&z->exp, &x->exp, bits - z->prec);
        }
    }
}
