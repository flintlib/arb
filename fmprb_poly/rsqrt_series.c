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
_fmprb_poly_rsqrt_series(fmprb_ptr g,
    fmprb_srcptr h, long hlen, long len, long prec)
{
    hlen = FLINT_MIN(hlen, len);

    fmprb_rsqrt(g, h, prec);

    if (hlen == 1)
    {
        _fmprb_vec_zero(g + 1, len - 1);
    }
    else if (len == 2)
    {
        fmprb_div(g + 1, h + 1, h, prec);
        fmprb_mul(g + 1, g + 1, g, prec);
        fmprb_mul_2exp_si(g + 1, g + 1, -1);
        fmprb_neg(g + 1, g + 1);
    }
    else
    {
        fmprb_ptr t, u;
        long tlen;
        t = _fmprb_vec_init(2 * len);
        u = t + len;

        NEWTON_INIT(1, len)

        NEWTON_LOOP(m, n)
        tlen = FLINT_MIN(2 * m - 1, n);
        _fmprb_poly_mullow(t, g, m, g, m, tlen, prec);
        _fmprb_poly_mullow(u, g, m, t, tlen, n, prec);
        _fmprb_poly_mullow(t, u, n, h, hlen, n, prec);
        _fmprb_vec_scalar_mul_2exp_si(g + m, t + m, n - m, -1);
        _fmprb_vec_neg(g + m, g + m, n - m);
        NEWTON_END_LOOP

        NEWTON_END

        _fmprb_vec_clear(t, 2 * len);
    }
}

void
fmprb_poly_rsqrt_series(fmprb_poly_t g, const fmprb_poly_t h, long n, long prec)
{
    if (n == 0)
    {
        fmprb_poly_zero(g);
        return;
    }

    if (g == h)
    {
        fmprb_poly_t t;
        fmprb_poly_init(t);
        fmprb_poly_rsqrt_series(t, h, n, prec);
        fmprb_poly_swap(g, t);
        fmprb_poly_clear(t);
        return;
    }

    fmprb_poly_fit_length(g, n);
    if (h->length == 0)
        _fmprb_vec_indeterminate(g->coeffs, n);
    else
        _fmprb_poly_rsqrt_series(g->coeffs, h->coeffs, h->length, n, prec);
    _fmprb_poly_set_length(g, n);
    _fmprb_poly_normalise(g);
}

