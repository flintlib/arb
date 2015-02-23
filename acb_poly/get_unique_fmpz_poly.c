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

#include "acb_poly.h"

int
acb_poly_get_unique_fmpz_poly(fmpz_poly_t res, const acb_poly_t src)
{
    int success;

    fmpz_poly_fit_length(res, src->length);
    success = _acb_vec_get_unique_fmpz_vec(res->coeffs, src->coeffs, src->length);
    _fmpz_poly_set_length(res, src->length);
    _fmpz_poly_normalise(res);
    return success;
}

