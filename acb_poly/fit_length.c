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
acb_poly_fit_length(acb_poly_t poly, long len)
{
    long i;

    if (len > poly->alloc)
    {
        if (len < 2 * poly->alloc)
            len = 2 * poly->alloc;

        poly->coeffs = flint_realloc(poly->coeffs,
            len * sizeof(acb_struct));

        for (i = poly->alloc; i < len; i++)
            acb_init(poly->coeffs + i);

        poly->alloc = len;
    }
}
