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

#include "zeta.h"
#include "hypgeom.h"

static void
atanh_bsplit(fmprb_t s, ulong p, ulong q, long prec)
{
    fmprb_t t;
    hypgeom_t series;
    hypgeom_init(series);
    fmprb_init(t);

    fmpz_poly_set_ui(series->A, 1);
    fmpz_poly_set_str(series->B, "2  1 2");
    fmpz_poly_set_ui(series->P, p);
    fmpz_mul(series->P->coeffs, series->P->coeffs, series->P->coeffs);
    fmpz_poly_set_ui(series->Q, q);
    fmpz_mul(series->Q->coeffs, series->Q->coeffs, series->Q->coeffs);

    fmprb_hypgeom_infsum(s, t, series, prec, prec);
    fmprb_mul_ui(s, s, p, prec);
    fmprb_mul_ui(t, t, q, prec);
    fmprb_div(s, s, t, prec);

    fmprb_clear(t);
    hypgeom_clear(series);
}

void
zeta_log_ui_from_prev(fmprb_t s, ulong k, fmprb_t log_prev, ulong prev, long prec)
{
    if (k > 200 && prec > 5000)
    {
        fmprb_t t;
        fmprb_init(t);

        atanh_bsplit(t, k - prev, k + prev, prec);
        fmprb_mul_2exp_si(t, t, 1);
        fmprb_add(s, log_prev, t, prec);

        fmprb_clear(t);
    }
    else
    {
        fmprb_log_ui(s, k, prec);
    }
}

