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
arb_poly_set_fmpq_poly(arb_poly_t poly, const fmpq_poly_t src)
{
    if (fmpq_poly_is_zero(src))
    {
        arb_poly_zero(poly);
    }
    else
    {
        long prec, len, numbits, denbits, shift;

        len = src->length;
        prec = poly->prec;

        _arb_poly_fit_length(poly, len);

        /* todo: handle power of two denominators exactly */

        numbits = _fmpz_vec_max_bits(fmpq_poly_numref(src), len);
        numbits = FLINT_ABS(numbits);
        denbits = fmpz_bits(fmpq_poly_denref(src));

        shift = prec - numbits + denbits;

        if (shift >= 0)
        {
            _fmpz_vec_scalar_mul_2exp(poly->coeffs,
                fmpq_poly_numref(src), len, shift);
            _fmpz_vec_scalar_tdiv_q_fmpz(poly->coeffs,
                poly->coeffs, len, fmpq_poly_denref(src));
        }
        else
        {
            fmpz_t d;
            fmpz_init(d);
            fmpz_mul_2exp(d, fmpq_poly_denref(src), -shift);
            _fmpz_vec_scalar_tdiv_q_fmpz(poly->coeffs, src->coeffs, len, d);
            fmpz_clear(d);
        }

        fmpz_one(&poly->rad);
        fmpz_set_si(&poly->exp, -shift);

        poly->length = len;
    }
}
