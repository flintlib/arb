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

#include "arb_poly.h"

void
_arb_poly_evaluate_acb(acb_t res, arb_srcptr f, long len,
                           const acb_t x, long prec)
{
    _arb_poly_evaluate_acb_rectangular(res, f, len, x, prec);
}

void
arb_poly_evaluate_acb(acb_t res, const arb_poly_t f, const acb_t a, long prec)
{
    _arb_poly_evaluate_acb(res, f->coeffs, f->length, a, prec);
}

