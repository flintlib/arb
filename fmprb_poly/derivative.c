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

#include "fmprb_poly.h"

void
_fmprb_poly_derivative(fmprb_struct * res, const fmprb_struct * poly, long len, long prec)
{
    long i;

    for (i = 1; i < len; i++)
        fmprb_mul_ui(res + i - 1, poly + i, i, prec);
}

void
fmprb_poly_derivative(fmprb_poly_t res, const fmprb_poly_t poly, long prec)
{
    long len = poly->length;

    if (len < 2)
    {
        fmprb_poly_zero(res);
    }
    else
    {
        fmprb_poly_fit_length(res, len - 1);
        _fmprb_poly_derivative(res->coeffs, poly->coeffs, len, prec);
        _fmprb_poly_set_length(res, len - 1);
    }
}