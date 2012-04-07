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

int
arb_poly_contains_fmpq_poly(const arb_poly_t x, const fmpq_poly_t y)
{
    long i, exp, xlen, ylen;

    xlen = x->length;
    ylen = y->length;

    if (xlen < ylen)
    {
        return 0;
    }

    exp = *arb_poly_expref(x);

    if (COEFF_IS_MPZ(exp))
    {
        printf("arb_poly_contains_fmpq_poly: huge exponent\n");
        abort();
    }

    for (i = 0; i < ylen; i++)
    {
        if (!_arb_contains_fmpq(arb_poly_coeffs(x) + i,
            arb_poly_radref(x), exp, fmpq_poly_numref(y) + i,
            fmpq_poly_denref(y)))
        {
            return 0;
        }
    }

    /* remaining intervals must contain zero */
    for (i = ylen; i < xlen; i++)
        if (!(fmpz_cmpabs(arb_poly_radref(x), arb_poly_coeffs(x) + i) >= 0))
            return 0;

    return 1;
}
