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

void
_arb_poly_add(arb_ptr res, arb_srcptr poly1, long len1,
    arb_srcptr poly2, long len2, long prec)
{
    long i, min = FLINT_MIN(len1, len2);

    for (i = 0; i < min; i++)
        arb_add(res + i, poly1 + i, poly2 + i, prec);

    for (i = min; i < len1; i++)
        arb_set_round(res + i, poly1 + i, prec);

    for (i = min; i < len2; i++)
        arb_set_round(res + i, poly2 + i, prec);
}

void
arb_poly_add(arb_poly_t res, const arb_poly_t poly1,
              const arb_poly_t poly2, long prec)
{
    long max = FLINT_MAX(poly1->length, poly2->length);

    arb_poly_fit_length(res, max);

    _arb_poly_add(res->coeffs, poly1->coeffs, poly1->length, poly2->coeffs,
                   poly2->length, prec);

    _arb_poly_set_length(res, max);
    _arb_poly_normalise(res);
}

