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
_arb_poly_set_coeff_same_exp(arb_poly_t z, long i, const arb_t c)
{
    long e1, e2, j;
    fmpz_t rad;

    e1 = z->exp;
    e2 = c->exp;

    if (z->length < i + 1)
    {
        _arb_poly_fit_length(z, i + 1);
        for (j = z->length; j < i; j++)
            fmpz_zero(z->coeffs + j);
        z->length = i + 1;
    }

    fmpz_init(rad);

    if (e1 > e2)
    {
        fmpz_tdiv_q_2exp(z->coeffs + i, arb_midref(c), e1 - e2);
        fmpz_cdiv_q_2exp(rad, arb_radref(c), e1 - e2);
        fmpz_add_ui(rad, rad, 1);
    }
    else
    {
        fmpz_mul_2exp(z->coeffs + i, arb_midref(c), e2 - e1);
        fmpz_mul_2exp(rad, arb_radref(c), e2 - e1);
    }

    if (fmpz_cmp(arb_poly_radref(z), rad) < 0)
        fmpz_set(arb_poly_radref(z), rad);

    fmpz_clear(rad);
}
