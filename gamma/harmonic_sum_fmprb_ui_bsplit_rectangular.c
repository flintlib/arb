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
bsplit_step(fmprb_t p, fmprb_t q, fmprb_srcptr poly, ulong step,
        const fmprb_t x, ulong a, ulong b, long prec)
{
    if (b - a == step)
    {
        fmprb_add_ui(p, x, a, prec);
        _fmprb_poly_evaluate2_rectangular(q, p, poly, step + 1, p, prec);
    }
    else
    {
        fmprb_t r, s;
        ulong m;

        fmprb_init(r);
        fmprb_init(s);

        m = a + ((b - a) / (2 * step)) * step;
        bsplit_step(p, q, poly, step, x, a, m, prec);
        bsplit_step(r, s, poly, step, x, m, b, prec);

        fmprb_mul(p, p, s, prec);
        fmprb_mul(r, r, q, prec);
        fmprb_add(p, p, r, prec);
        fmprb_mul(q, q, s, prec);

        fmprb_clear(r);
        fmprb_clear(s);
    }
}

void
gamma_harmonic_sum_fmprb_ui_bsplit_rectangular(fmprb_t y, const fmprb_t x, ulong n, ulong step, long prec)
{
    fmprb_t t, u;
    ulong b, s1, s2;
    long wp;

    wp = FMPR_PREC_ADD(prec, FLINT_BIT_COUNT(n));

    fmprb_init(t);
    fmprb_init(u);

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
        _fmprb_poly_rising_ui_series(poly, xpoly, 2, step, step + 1, wp);
        bsplit_step(t, u, poly, step, x, 0, b, wp);
        fmprb_div(t, t, u, wp);
        _fmprb_vec_clear(poly, step + 1);
    }
    else
    {
        fmprb_zero(t);
    }

    fmprb_add_ui(u, x, b, wp);
    gamma_harmonic_sum_fmprb_ui_bsplit(u, u, n - b, wp);
    fmprb_add(y, t, u, prec);

    fmprb_clear(t);
    fmprb_clear(u);
}

