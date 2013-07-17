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

#define SIN_TAN_CUTOFF 20

void
_fmprb_poly_sin_series(fmprb_ptr g, const fmprb_srcptr h, long hlen, long len, long prec)
{
    fmprb_ptr t, u;
    hlen = FLINT_MIN(hlen, len);

    if (hlen == 1)
    {
        fmprb_sin(g, h, prec);
        _fmprb_vec_zero(g + 1, len - 1);
    }
    else if (hlen < SIN_TAN_CUTOFF)
    {
        t = _fmprb_vec_init(len);
        _fmprb_poly_sin_cos_series_basecase(g, t, h, hlen, len, prec);
        _fmprb_vec_clear(t, len);
    }
    else
    {
        /*
            sin(x) = 2*tan(x/2)/(1+tan(x/2)^2)
            cos(x) = (1-tan(x/2)^2)/(1+tan(x/2)^2)

            In order to use the tangent expansion at the origin:
            sin(h0 + h1) = cos(h0) sin(h1) + sin(h0) cos(h1)
        */
        fmprb_ptr v;

        t = _fmprb_vec_init(3 * len);
        u = t + len;
        v = u + len;

        /* t = tan(h/2)^2 */
        _fmprb_vec_scalar_mul_2exp_si(t, h, hlen, -1);
        fmprb_zero(t);
        _fmprb_poly_tan_series(u, t, hlen, len, prec);
        _fmprb_poly_mullow(t, u, len, u, len, len, prec);

        /* v = 1/(1 + t) */
        fmprb_add_ui(t, t, 1, prec);
        _fmprb_poly_inv_series(v, t, len, len, prec);



        _fmprb_vec_scalar_mul_2exp_si(g, g, len, 1);
        _fmprb_vec_clear(t, 2 * len);
    }
}

void
fmprb_poly_sin_series(fmprb_poly_t g, const fmprb_poly_t h, long n, long prec)
{
    if (h->length == 0 || n == 0)
    {
        fmprb_poly_zero(g);
        return;
    }

    fmprb_poly_fit_length(g, n);
    _fmprb_poly_sin_series(g->coeffs, h->coeffs, h->length, n, prec);
    _fmprb_poly_set_length(g, n);
    _fmprb_poly_normalise(g);
}

