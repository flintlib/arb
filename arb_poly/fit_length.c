/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "arb_poly.h"

void
_arb_poly_fit_length(arb_poly_t poly, long length)
{
    if (length > poly->alloc)
    {
        long i, alloc = FLINT_MAX(length, 2 * poly->alloc);

        poly->coeffs = (fmpz *) flint_realloc(poly->coeffs,
            alloc * sizeof(fmpz));
        for (i = poly->alloc; i < alloc; i++)
            poly->coeffs[i] = 0;

        poly->alloc = alloc;
    }
}
