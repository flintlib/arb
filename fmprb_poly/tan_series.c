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

#define TAN_NEWTON_CUTOFF 20

void
_fmprb_poly_tan_series(fmprb_struct * g,
    const fmprb_struct * h, long hlen, long len, long prec)
{
    fmprb_struct *t, *u;
    hlen = FLINT_MIN(hlen, len);

    if (hlen == 1)
    {
        fmprb_tan(g, h, prec);
        _fmprb_vec_zero(g + 1, len - 1);
        return;
    }

    t = _fmprb_vec_init(2 * len);
    u = t + len;

    NEWTON_INIT(TAN_NEWTON_CUTOFF, len)

    NEWTON_BASECASE(n)
    _fmprb_poly_sin_cos_series_basecase(t, u, h, hlen, n, prec);
    _fmprb_poly_div_series(g, t, len, u, n, n, prec);
    NEWTON_END_BASECASE

    NEWTON_LOOP(m, n)
    _fmprb_poly_mullow(u, g, m, g, m, n, prec);
    fmprb_add_ui(u, u, 1, prec);
    _fmprb_poly_atan_series(t, g, m, n, prec);
    _fmprb_poly_sub(t + m, h + m, FLINT_MAX(0, hlen - m), t + m, n - m, prec);
    _fmprb_poly_mullow(g + m, u, n, t + m, n - m, n - m, prec);
    NEWTON_END_LOOP

    NEWTON_END

    _fmprb_vec_clear(t, 2 * len);
}

void
fmprb_poly_tan_series(fmprb_poly_t g, const fmprb_poly_t h, long n, long prec)
{
    if (h->length == 0 || n == 0)
    {
        fmprb_poly_zero(g);
        return;
    }

    if (g == h)
    {
        fmprb_poly_t t;
        fmprb_poly_init(t);
        fmprb_poly_tan_series(t, h, n, prec);
        fmprb_poly_swap(g, t);
        fmprb_poly_clear(t);
        return;
    }

    fmprb_poly_fit_length(g, n);
    _fmprb_poly_tan_series(g->coeffs, h->coeffs, h->length, n, prec);
    _fmprb_poly_set_length(g, n);
    _fmprb_poly_normalise(g);
}

