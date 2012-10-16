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

#include "fmpcb_poly.h"

void
fmpcb_poly_set_fmpq_poly(fmpcb_poly_t poly, const fmpq_poly_t re, long prec)
{
    fmprb_poly_t t;
    fmprb_poly_init(t);
    fmprb_poly_set_fmpq_poly(t, re, prec);
    fmpcb_poly_set_fmprb_poly(poly, t);
    fmprb_poly_clear(t);
}

void
fmpcb_poly_set2_fmpq_poly(fmpcb_poly_t poly, const fmpq_poly_t re, const fmpq_poly_t im, long prec)
{
    fmprb_poly_t t, u;
    fmprb_poly_init(t);
    fmprb_poly_init(u);

    fmprb_poly_set_fmpq_poly(t, re, prec);
    fmprb_poly_set_fmpq_poly(u, im, prec);

    fmpcb_poly_set2_fmprb_poly(poly, t, u);

    fmprb_poly_clear(t);
    fmprb_poly_clear(u);
}
