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

#include "acb_poly.h"

void
_acb_poly_derivative(acb_ptr res, acb_srcptr poly, long len, long prec)
{
    long i;

    for (i = 1; i < len; i++)
        acb_mul_ui(res + i - 1, poly + i, i, prec);
}

void
acb_poly_derivative(acb_poly_t res, const acb_poly_t poly, long prec)
{
    long len = poly->length;

    if (len < 2)
    {
        acb_poly_zero(res);
    }
    else
    {
        acb_poly_fit_length(res, len - 1);
        _acb_poly_derivative(res->coeffs, poly->coeffs, len, prec);
        _acb_poly_set_length(res, len - 1);
    }
}
