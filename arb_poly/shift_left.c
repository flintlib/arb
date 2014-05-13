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

    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include "arb_poly.h"

void
_arb_poly_shift_left(arb_ptr res, arb_srcptr poly, long len, long n)
{
    long i;

    /* Copy in reverse to avoid writing over unshifted coefficients */
    if (res != poly)
    {
        for (i = len; i--; )
            arb_set(res + n + i, poly + i);
    }
    else
    {
        for (i = len; i--; )
            arb_swap(res + n + i, res + i);
    }

    for (i = 0; i < n; i++)
        arb_zero(res + i);
}

void
arb_poly_shift_left(arb_poly_t res, const arb_poly_t poly, long n)
{
    if (n == 0)
    {
        arb_poly_set(res, poly);
        return;
    }

    if (poly->length == 0)
    {
        arb_poly_zero(res);
        return;
    }

    arb_poly_fit_length(res, poly->length + n);
    _arb_poly_shift_left(res->coeffs, poly->coeffs, poly->length, n);
    _arb_poly_set_length(res, poly->length + n);
}

