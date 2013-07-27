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


void
_fmpcb_poly_atan_series(fmpcb_ptr g, fmpcb_srcptr h, long hlen, long n, long prec)
{
    fmpcb_t c;
    fmpcb_init(c);

    /* TODO: fmpcb_atan(c, h, prec); */
    if (!fmpcb_is_zero(h))
        abort();
    fmpcb_zero(g);

    hlen = FLINT_MIN(hlen, n);

    if (hlen == 1)
    {
        _fmpcb_vec_zero(g + 1, n - 1);
    }
    else
    {
        fmpcb_ptr t, u;
        long ulen;

        t = _fmpcb_vec_init(n);
        u = _fmpcb_vec_init(n);

        /* atan(h(x)) = integral(h'(x)/(1+h(x)^2)) */
        ulen = FLINT_MIN(n, 2 * hlen - 1);
        _fmpcb_poly_mullow(u, h, hlen, h, hlen, ulen, prec);
        fmpcb_add_ui(u, u, 1, prec);

        _fmpcb_poly_derivative(t, h, hlen, prec);
        _fmpcb_poly_div_series(g, t, hlen - 1, u, ulen, n, prec);
        _fmpcb_poly_integral(g, g, n, prec);

        _fmpcb_vec_clear(t, n);
        _fmpcb_vec_clear(u, n);
    }

    fmpcb_swap(g, c);
    fmpcb_clear(c);
}

void
fmpcb_poly_atan_series(fmpcb_poly_t g, const fmpcb_poly_t h, long n, long prec)
{
    long hlen = h->length;

    if (hlen == 0 || n == 0)
    {
        fmpcb_poly_zero(g);
        return;
    }

    fmpcb_poly_fit_length(g, n);
    _fmpcb_poly_atan_series(g->coeffs, h->coeffs, hlen, n, prec);
    _fmpcb_poly_set_length(g, n);
    _fmpcb_poly_normalise(g);
}

