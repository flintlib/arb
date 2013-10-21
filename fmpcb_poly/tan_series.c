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

#include "fmpcb_poly.h"

#define TAN_NEWTON_CUTOFF 20

void
_fmpcb_poly_tan_series(fmpcb_ptr g,
    fmpcb_srcptr h, long hlen, long len, long prec)
{
    hlen = FLINT_MIN(hlen, len);

    if (hlen == 1)
    {
        fmpcb_tan(g, h, prec);
        _fmpcb_vec_zero(g + 1, len - 1);
    }
    else if (len == 2)
    {
        fmpcb_t t;
        fmpcb_init(t);
        fmpcb_tan(g, h, prec);
        fmpcb_mul(t, g, g, prec);
        fmpcb_add_ui(t, t, 1, prec);
        fmpcb_mul(g + 1, t, h + 1, prec);  /* safe since hlen >= 2 */
        fmpcb_clear(t);
    }
    else
    {
        fmpcb_ptr t, u;

        t = _fmpcb_vec_init(2 * len);
        u = t + len;

        NEWTON_INIT(TAN_NEWTON_CUTOFF, len)

        NEWTON_BASECASE(n)
        _fmpcb_poly_sin_cos_series_basecase(t, u, h, hlen, n, prec);
        _fmpcb_poly_div_series(g, t, n, u, n, n, prec);
        NEWTON_END_BASECASE

        NEWTON_LOOP(m, n)
        _fmpcb_poly_mullow(u, g, m, g, m, n, prec);
        fmpcb_add_ui(u, u, 1, prec);
        _fmpcb_poly_atan_series(t, g, m, n, prec);
        _fmpcb_poly_sub(t + m, h + m, FLINT_MAX(0, hlen - m), t + m, n - m, prec);
        _fmpcb_poly_mullow(g + m, u, n, t + m, n - m, n - m, prec);
        NEWTON_END_LOOP

        NEWTON_END

        _fmpcb_vec_clear(t, 2 * len);
    }
}

void
fmpcb_poly_tan_series(fmpcb_poly_t g, const fmpcb_poly_t h, long n, long prec)
{
    if (h->length == 0 || n == 0)
    {
        fmpcb_poly_zero(g);
        return;
    }

    if (g == h)
    {
        fmpcb_poly_t t;
        fmpcb_poly_init(t);
        fmpcb_poly_tan_series(t, h, n, prec);
        fmpcb_poly_swap(g, t);
        fmpcb_poly_clear(t);
        return;
    }

    fmpcb_poly_fit_length(g, n);
    _fmpcb_poly_tan_series(g->coeffs, h->coeffs, h->length, n, prec);
    _fmpcb_poly_set_length(g, n);
    _fmpcb_poly_normalise(g);
}

