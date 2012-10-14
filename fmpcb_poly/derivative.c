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

void
_fmpcb_poly_derivative(fmpcb_struct * res, const fmpcb_struct * poly, long len, long prec)
{
    long i;

    for (i = 1; i < len; i++)
        fmpcb_mul_ui(res + i - 1, poly + i, i, prec);
}

void
fmpcb_poly_derivative(fmpcb_poly_t res, const fmpcb_poly_t poly, long prec)
{
    long len = poly->length;

    if (len < 2)
    {
        fmpcb_poly_zero(res);
    }
    else
    {
        fmpcb_poly_fit_length(res, len - 1);
        _fmpcb_poly_derivative(res->coeffs, poly->coeffs, len, prec);
        _fmpcb_poly_set_length(res, len - 1);
    }
}
