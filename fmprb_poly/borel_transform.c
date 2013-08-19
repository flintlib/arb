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

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include "fmprb_poly.h"

void
_fmprb_poly_borel_transform(fmprb_ptr res, fmprb_srcptr poly, long len, long prec)
{
    long i;

    fmprb_t t;
    fmprb_init(t);

    fmprb_one(t);

    for (i = 0; i < len; i++)
    {
        if (i > 1)
            fmprb_mul_ui(t, t, i, prec);

        fmprb_div(res + i, poly + i, t, prec);
    }

    fmprb_clear(t);
}

void
fmprb_poly_borel_transform(fmprb_poly_t res, const fmprb_poly_t poly, long prec)
{
    fmprb_poly_fit_length(res, poly->length);
    _fmprb_poly_borel_transform(res->coeffs, poly->coeffs, poly->length, prec);
    _fmprb_poly_set_length(res, poly->length);
    _fmprb_poly_normalise(res);
}

