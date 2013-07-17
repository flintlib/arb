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
#include "fmprb_poly.h"

static void
bsplit_step(fmpcb_t y, fmprb_srcptr poly, ulong step, const fmpcb_t x, ulong a, ulong b, long prec)
{
    if (b - a == step)
    {
        fmpcb_add_ui(y, x, a, prec);
        _fmprb_poly_evaluate_fmpcb(y, poly, step + 1, y, prec);
    }
    else
    {
        ulong m = a + ((b - a) / (2 * step)) * step;
        fmpcb_t t;
        fmpcb_init(t);
        bsplit_step(y, poly, step, x, a, m, prec);
        bsplit_step(t, poly, step, x, m, b, prec);
        fmpcb_mul(y, y, t, prec);
        fmpcb_clear(t);
    }
}

void
gamma_rising_fmpcb_ui_bsplit_rectangular(fmpcb_t y, const fmpcb_t x, ulong n, ulong step, long prec)
{
    fmpcb_t t, u;
    ulong s1, s2, b;
    long wp;

    wp = FMPR_PREC_ADD(prec, FLINT_BIT_COUNT(n));

    fmpcb_init(t);
    fmpcb_init(u);

    if (step == 0)
    {
        s1 = 2 * sqrt(n);
        s2 = 10 * pow(prec, 0.25);
        step = FLINT_MIN(s1, s2);
    }
    step = FLINT_MAX(step, 1);
    b = (n / step) * step;

    if (b != 0)
    {
        fmprb_ptr poly;
        fmprb_struct xpoly[2];
        poly = _fmprb_vec_init(step + 1);
        fmprb_init(xpoly + 0);
        fmprb_init(xpoly + 1);
        fmprb_one(xpoly + 1);
        _fmprb_poly_rfac_series_ui(poly, xpoly, 2, step, step + 1, wp);
        bsplit_step(t, poly, step, x, 0, b, wp);
        _fmprb_vec_clear(poly, step + 1);
    }
    else
    {
        fmpcb_one(t);
    }

    fmpcb_add_ui(u, x, b, wp);
    gamma_rising_fmpcb_ui_bsplit_simple(u, u, n - b, wp);
    fmpcb_mul(y, t, u, prec);

    fmpcb_clear(t);
    fmpcb_clear(u);
}

