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
_arb_poly_asin_series(arb_ptr g, arb_srcptr h, long hlen, long n, long prec)
{
    arb_t c;
    arb_init(c);

    arb_asin(c, h, prec);

    hlen = FLINT_MIN(hlen, n);

    if (hlen == 1)
    {
        _arb_vec_zero(g + 1, n - 1);
    }
    else
    {
        arb_ptr t, u;
        long ulen;

        t = _arb_vec_init(n);
        u = _arb_vec_init(n);

        /* asin(h(x)) = integral(h'(x)/sqrt(1-h(x)^2)) */
        ulen = FLINT_MIN(n, 2 * hlen - 1);
        _arb_poly_mullow(u, h, hlen, h, hlen, ulen, prec);
        arb_sub_ui(u, u, 1, prec);
        _arb_vec_neg(u, u, ulen);
        _arb_poly_rsqrt_series(t, u, ulen, n, prec);
        _arb_poly_derivative(u, h, hlen, prec);
        _arb_poly_mullow(g, t, n, u, hlen - 1, n, prec);
        _arb_poly_integral(g, g, n, prec);

        _arb_vec_clear(t, n);
        _arb_vec_clear(u, n);
    }

    arb_swap(g, c);
    arb_clear(c);
}

void
arb_poly_asin_series(arb_poly_t g, const arb_poly_t h, long n, long prec)
{
    long hlen = h->length;

    if (hlen == 0 || n == 0)
    {
        arb_poly_zero(g);
        return;
    }

    arb_poly_fit_length(g, n);
    _arb_poly_asin_series(g->coeffs, h->coeffs, hlen, n, prec);
    _arb_poly_set_length(g, n);
    _arb_poly_normalise(g);
}

