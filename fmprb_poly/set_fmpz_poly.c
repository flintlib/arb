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

/* XXX: refactor this */
void fmprb_set_fmpz_round(fmprb_t y, const fmpz_t x, long prec)
{
    fmprb_set_fmpz(y, x);

    if (!fmpz_is_zero(x))
    {
        long r;

        r = fmpr_set_round(fmprb_midref(y), fmprb_midref(y), prec, FMPR_RND_DOWN);
        fmpr_set_error_result(fmprb_radref(y), fmprb_midref(y), r);
    }
}



void
fmprb_poly_set_fmpz_poly(fmprb_poly_t poly, const fmpz_poly_t src, long prec)
{
    long i, len = fmpz_poly_length(src);

    fmprb_poly_fit_length(poly, len);
    _fmprb_poly_set_length(poly, len);

    for (i = 0; i < len; i++)
        fmprb_set_fmpz_round(poly->coeffs + i, src->coeffs + i, prec);
}
