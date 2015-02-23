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

    Copyright (C) 2015 Fredrik Johansson

******************************************************************************/

#include "arb_poly.h"

void
arb_poly_add_si(arb_poly_t res, const arb_poly_t x, long y, long prec)
{
    long len = x->length;

    if (len == 0)
    {
        arb_poly_set_si(res, y);
    }
    else
    {
        arb_poly_fit_length(res, len);

        arb_add_si(res->coeffs, x->coeffs, y, prec);

        if (res != x)
            _arb_vec_set(res->coeffs + 1, x->coeffs + 1, len - 1);

        _arb_poly_set_length(res, len);
        _arb_poly_normalise(res);
    }
}

