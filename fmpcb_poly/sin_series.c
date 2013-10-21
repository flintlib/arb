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
_fmpcb_poly_sin_series(fmpcb_ptr g, fmpcb_srcptr h, long hlen, long n, long prec)
{
    hlen = FLINT_MIN(hlen, n);

    if (hlen == 1)
    {
        fmpcb_sin(g, h, prec);
        _fmpcb_vec_zero(g + 1, n - 1);
    }
    else if (n == 2)
    {
        fmpcb_t t;
        fmpcb_init(t);
        fmpcb_sin_cos(g, t, h, prec);
        fmpcb_mul(g + 1, h + 1, t, prec);  /* safe since hlen >= 2 */
        fmpcb_clear(t);
    }
    else
    {
        fmpcb_ptr t = _fmpcb_vec_init(n);
        _fmpcb_poly_sin_cos_series(g, t, h, hlen, n, prec);
        _fmpcb_vec_clear(t, n);
    }
}

void
fmpcb_poly_sin_series(fmpcb_poly_t g, const fmpcb_poly_t h, long n, long prec)
{
    long hlen = h->length;

    if (hlen == 0 || n == 0)
    {
        fmpcb_poly_zero(g);
        return;
    }

    if (hlen == 1)
        n = 1;

    fmpcb_poly_fit_length(g, n);
    _fmpcb_poly_sin_series(g->coeffs, h->coeffs, hlen, n, prec);
    _fmpcb_poly_set_length(g, n);
    _fmpcb_poly_normalise(g);
}

