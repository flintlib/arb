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

#include "arb.h"
#include "hypgeom.h"

static void
atanh_bsplit(arb_t s, ulong c, long a, long prec)
{
    arb_t t;
    hypgeom_t series;
    hypgeom_init(series);
    arb_init(t);

    fmpz_poly_set_ui(series->A, 1);
    fmpz_poly_set_coeff_ui(series->B, 0, 1);
    fmpz_poly_set_coeff_ui(series->B, 1, 2);
    fmpz_poly_set_ui(series->P, 1);
    fmpz_poly_set_ui(series->Q, c * c);

    arb_hypgeom_infsum(s, t, series, prec, prec);
    arb_mul_si(s, s, a, prec);
    arb_mul_ui(t, t, c, prec);
    arb_div(s, s, t, prec);

    arb_clear(t);
    hypgeom_clear(series);
}

void
arb_const_log10_eval(arb_t s, long prec)
{
    arb_t t;
    arb_init(t);

    prec += FLINT_CLOG2(prec);

    atanh_bsplit(s, 31, 46, prec);
    atanh_bsplit(t, 49, 34, prec);
    arb_add(s, s, t, prec);
    atanh_bsplit(t, 161, 20, prec);
    arb_add(s, s, t, prec);

    arb_clear(t);
}

ARB_DEF_CACHED_CONSTANT(arb_const_log10, arb_const_log10_eval)

