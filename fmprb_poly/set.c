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
fmprb_poly_set(fmprb_poly_t poly, const fmprb_poly_t src)
{
    long i, len = fmprb_poly_length(src);

    fmprb_poly_fit_length(poly, len);
    _fmprb_poly_set_length(poly, len);

    for (i = 0; i < len; i++)
        fmprb_set(poly->coeffs + i, src->coeffs + i);
}
