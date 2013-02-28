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

#include "gamma.h"
#include "hypgeom.h"

void
gamma_const_1_3_eval(fmprb_t s, long prec)
{
    hypgeom_t series;
    fmprb_t t, u;

    fmprb_init(t);
    fmprb_init(u);

    hypgeom_init(series);

    fmpz_poly_set_str(series->A, "1  1");
    fmpz_poly_set_str(series->B, "1  1");
    fmpz_poly_set_str(series->P, "4  5 -46 108 -72");
    fmpz_poly_set_str(series->Q, "4  0 0 0 512000");

    prec += FLINT_CLOG2(prec);
    fmprb_hypgeom_infsum(s, t, series, prec, prec);

    fmprb_sqrt_ui(u, 10, prec);
    fmprb_mul(t, t, u, prec);

    fmprb_const_pi(u, prec);
    fmprb_pow_ui(u, u, 4, prec);
    fmprb_mul_ui(u, u, 12, prec);
    fmprb_mul(s, s, u, prec);

    fmprb_div(s, s, t, prec);
    fmprb_root(s, s, 2, prec);
    fmprb_root(s, s, 3, prec);

    hypgeom_clear(series);
    fmprb_clear(t);
    fmprb_clear(u);
}

DEF_CACHED_CONSTANT(gamma_const_1_3, gamma_const_1_3_eval)

