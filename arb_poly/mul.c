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

#include "arb_poly.h"

void _arb_poly_mul(arb_ptr C,
    arb_srcptr A, long lenA,
    arb_srcptr B, long lenB, long prec)
{
    _arb_poly_mullow(C, A, lenA, B, lenB, lenA + lenB - 1, prec);
}

void
arb_poly_mul(arb_poly_t res, const arb_poly_t poly1,
              const arb_poly_t poly2, long prec)
{
    long len_out;

    if ((poly1->length == 0) || (poly2->length == 0))
    {
        arb_poly_zero(res);
        return;
    }

    len_out = poly1->length + poly2->length - 1;

    if (res == poly1 || res == poly2)
    {
        arb_poly_t temp;
        arb_poly_init2(temp, len_out);
        _arb_poly_mul(temp->coeffs, poly1->coeffs, poly1->length,
                                 poly2->coeffs, poly2->length, prec);
        arb_poly_swap(res, temp);
        arb_poly_clear(temp);
    }
    else
    {
        arb_poly_fit_length(res, len_out);
        _arb_poly_mul(res->coeffs, poly1->coeffs, poly1->length,
                                 poly2->coeffs, poly2->length, prec);
    }

    _arb_poly_set_length(res, len_out);
    _arb_poly_normalise(res);
}
