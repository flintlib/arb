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
_fmprb_poly_evaluate_fmpcb(fmpcb_t res, const fmprb_struct * f, long len,
                           const fmpcb_t x, long prec)
{
    _fmprb_poly_evaluate_fmpcb_rectangular(res, f, len, x, prec);
}

void
fmprb_poly_evaluate_fmpcb(fmpcb_t res, const fmprb_poly_t f, const fmpcb_t a, long prec)
{
    _fmprb_poly_evaluate_fmpcb(res, f->coeffs, f->length, a, prec);
}

