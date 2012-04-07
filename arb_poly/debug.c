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
arb_poly_debug(const arb_poly_t x)
{
    long i;

    printf("arb_poly(length=%ld, coeffs=[", x->length);
    for (i = 0; i < x->length; i++)
    {
        fmpz_print(arb_poly_coeffs(x) + i);
        if (i < x->length - 1)
            printf(", ");
    }
    printf("], rad=");
    fmpz_print(arb_poly_radref(x));
    printf(", exp=");
    fmpz_print(arb_poly_expref(x));
    printf(", prec=");
    printf("%ld", arb_poly_prec(x));
    printf(")");
}
