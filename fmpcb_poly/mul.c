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

#include "fmpcb_poly.h"

void _fmpcb_poly_mul(fmpcb_struct * C,
    const fmpcb_struct * A, long lenA,
    const fmpcb_struct * B, long lenB, long prec)
{
    _fmpcb_poly_mullow(C, A, lenA, B, lenB, lenA + lenB - 1, prec);
}

void
fmpcb_poly_mul(fmpcb_poly_t res, const fmpcb_poly_t poly1,
              const fmpcb_poly_t poly2, long prec)
{
    long len_out;

    if ((poly1->length == 0) || (poly2->length == 0))
    {
        fmpcb_poly_zero(res);
        return;
    }

    len_out = poly1->length + poly2->length - 1;

    if (res == poly1 || res == poly2)
    {
        fmpcb_poly_t temp;
        fmpcb_poly_init2(temp, len_out);
        _fmpcb_poly_mul(temp->coeffs, poly1->coeffs, poly1->length,
                                 poly2->coeffs, poly2->length, prec);
        fmpcb_poly_swap(res, temp);
        fmpcb_poly_clear(temp);
    }
    else
    {
        fmpcb_poly_fit_length(res, len_out);
        _fmpcb_poly_mul(res->coeffs, poly1->coeffs, poly1->length,
                                 poly2->coeffs, poly2->length, prec);
    }

    _fmpcb_poly_set_length(res, len_out);
    _fmpcb_poly_normalise(res);
}

