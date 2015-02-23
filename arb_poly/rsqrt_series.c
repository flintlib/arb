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

#include "arb_poly.h"

void
_arb_poly_rsqrt_series(arb_ptr g,
    arb_srcptr h, long hlen, long len, long prec)
{
    hlen = FLINT_MIN(hlen, len);

    while (hlen > 0 && arb_is_zero(h + hlen - 1))
        hlen--;

    if (hlen <= 1)
    {
        arb_rsqrt(g, h, prec);
        _arb_vec_zero(g + 1, len - 1);
    }
    else if (len == 2)
    {
        arb_rsqrt(g, h, prec);
        arb_div(g + 1, h + 1, h, prec);
        arb_mul(g + 1, g + 1, g, prec);
        arb_mul_2exp_si(g + 1, g + 1, -1);
        arb_neg(g + 1, g + 1);
    }
    else if (_arb_vec_is_zero(h + 1, hlen - 2))
    {
        arb_t t;
        arb_init(t);
        arf_set_si_2exp_si(arb_midref(t), -1, -1);
        _arb_poly_binomial_pow_arb_series(g, h, hlen, t, len, prec);
        arb_clear(t);
    }
    else
    {
        arb_ptr t, u;
        long tlen;
        t = _arb_vec_init(2 * len);
        u = t + len;

        arb_rsqrt(g, h, prec);

        NEWTON_INIT(1, len)

        NEWTON_LOOP(m, n)
        tlen = FLINT_MIN(2 * m - 1, n);
        _arb_poly_mullow(t, g, m, g, m, tlen, prec);
        _arb_poly_mullow(u, g, m, t, tlen, n, prec);
        _arb_poly_mullow(t, u, n, h, hlen, n, prec);
        _arb_vec_scalar_mul_2exp_si(g + m, t + m, n - m, -1);
        _arb_vec_neg(g + m, g + m, n - m);
        NEWTON_END_LOOP

        NEWTON_END

        _arb_vec_clear(t, 2 * len);
    }
}

void
arb_poly_rsqrt_series(arb_poly_t g, const arb_poly_t h, long n, long prec)
{
    if (n == 0)
    {
        arb_poly_zero(g);
        return;
    }

    if (g == h)
    {
        arb_poly_t t;
        arb_poly_init(t);
        arb_poly_rsqrt_series(t, h, n, prec);
        arb_poly_swap(g, t);
        arb_poly_clear(t);
        return;
    }

    arb_poly_fit_length(g, n);
    if (h->length == 0)
        _arb_vec_indeterminate(g->coeffs, n);
    else
        _arb_poly_rsqrt_series(g->coeffs, h->coeffs, h->length, n, prec);
    _arb_poly_set_length(g, n);
    _arb_poly_normalise(g);
}

