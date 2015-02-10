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
_arb_poly_majorant(arb_ptr res, arb_srcptr vec, long len, long prec)
{
    long i;

    for (i = 0; i < len; i++)
    {
        arb_get_abs_ubound_arf(arb_midref(res + i), vec + i, prec);
        mag_zero(arb_radref(res + i));
    }
}

void
arb_poly_majorant(arb_poly_t res, const arb_poly_t poly, long prec)
{
    arb_poly_fit_length(res, poly->length);
    _arb_poly_majorant(res->coeffs, poly->coeffs, poly->length, prec);
    _arb_poly_set_length(res, poly->length);
}

