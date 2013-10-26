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
_fmprb_poly_cos_series(fmprb_ptr g, fmprb_srcptr h, long hlen, long n, long prec)
{
    hlen = FLINT_MIN(hlen, n);

    if (hlen == 1)
    {
        fmprb_cos(g, h, prec);
        _fmprb_vec_zero(g + 1, n - 1);
    }
    else if (n == 2)
    {
        fmprb_t t;
        fmprb_init(t);
        fmprb_sin_cos(t, g, h, prec);
        fmprb_neg(t, t);
        fmprb_mul(g + 1, h + 1, t, prec);  /* safe since hlen >= 2 */
        fmprb_clear(t);
    }
    else
    {
        fmprb_ptr t = _fmprb_vec_init(n);
        _fmprb_poly_sin_cos_series(t, g, h, hlen, n, prec);
        _fmprb_vec_clear(t, n);
    }
}

void
fmprb_poly_cos_series(fmprb_poly_t g, const fmprb_poly_t h, long n, long prec)
{
    long hlen = h->length;

    if (n == 0)
    {
        fmprb_poly_zero(g);
        return;
    }

    if (hlen == 0)
    {
        fmprb_poly_one(g);
        return;
    }

    if (hlen == 1)
        n = 1;

    fmprb_poly_fit_length(g, n);
    _fmprb_poly_cos_series(g->coeffs, h->coeffs, hlen, n, prec);
    _fmprb_poly_set_length(g, n);
    _fmprb_poly_normalise(g);
}

