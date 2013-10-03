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
_fmpcb_poly_sqrt_series(fmpcb_ptr g,
    fmpcb_srcptr h, long hlen, long len, long prec)
{
    fmpcb_ptr t;
    hlen = FLINT_MIN(hlen, len);
    t = _fmpcb_vec_init(len);
    _fmpcb_poly_rsqrt_series(t, h, hlen, len, prec);
    _fmpcb_poly_mullow(g, t, len, h, hlen, len, prec);
    _fmpcb_vec_clear(t, len);
}

void
fmpcb_poly_sqrt_series(fmpcb_poly_t g, const fmpcb_poly_t h, long n, long prec)
{
    if (n == 0 || h->length == 0)
    {
        printf("fmpcb_poly_sqrt_series: require n > 0 and nonzero input\n");
        abort();
    }

    if (g == h)
    {
        fmpcb_poly_t t;
        fmpcb_poly_init(t);
        fmpcb_poly_sqrt_series(t, h, n, prec);
        fmpcb_poly_swap(g, t);
        fmpcb_poly_clear(t);
        return;
    }

    fmpcb_poly_fit_length(g, n);
    _fmpcb_poly_sqrt_series(g->coeffs, h->coeffs, h->length, n, prec);
    _fmpcb_poly_set_length(g, n);
    _fmpcb_poly_normalise(g);
}

