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

#include "fmprb.h"
#include "hypgeom.h"

void
fmprb_const_pi_chudnovsky(fmprb_t s, long prec)
{
    hypgeom_t series;
    fmprb_t t, u;

    fmprb_init(t);
    fmprb_init(u);
    hypgeom_init(series);

    fmpz_poly_set_str(series->A, "2  13591409 545140134");
    fmpz_poly_set_str(series->B, "1  1");
    fmpz_poly_set_str(series->P, "4  5 -46 108 -72");
    fmpz_poly_set_str(series->Q, "4  0 0 0 10939058860032000");

    prec += FLINT_CLOG2(prec) + 5;
    fmprb_hypgeom_infsum(s, t, series, prec, prec);

    fmprb_rsqrt_ui(u, 640320, prec);
    fmprb_mul(s, s, u, prec);

    fmprb_mul_ui(t, t, 640320 / 12, prec);
    fmprb_div(s, t, s, prec);

    hypgeom_clear(series);
    fmprb_clear(t);
    fmprb_clear(u);
}
