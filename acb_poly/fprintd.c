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
    Copyright (C) 2015 Arb authors

******************************************************************************/

#include "acb_poly.h"

void
acb_poly_fprintd(FILE * file, const acb_poly_t poly, slong digits)
{
    slong i;

    flint_fprintf(file, "[");

    for (i = 0; i < poly->length; i++)
    {
        acb_fprintd(file, poly->coeffs + i, digits);
        if (i + 1 < poly->length)
            flint_fprintf(file, "\n");
    }

    flint_fprintf(file, "]");
}
