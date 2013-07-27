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
_fmpcb_poly_sin_cos_series_tangent(fmpcb_ptr s, fmpcb_ptr c,
                        const fmpcb_srcptr h, long hlen, long len, long prec)
{
    fmpcb_ptr t, u, v;
    fmpcb_t s0, c0;
    hlen = FLINT_MIN(hlen, len);

    if (hlen == 1)
    {
        fmpcb_sin_cos(s, c, h, prec);
        _fmpcb_vec_zero(s + 1, len - 1);
        _fmpcb_vec_zero(c + 1, len - 1);
        return;
    }

    /*
    sin(x) = 2*tan(x/2)/(1+tan(x/2)^2)
    cos(x) = (1-tan(x/2)^2)/(1+tan(x/2)^2)
    */

    fmpcb_init(s0);
    fmpcb_init(c0);

    t = _fmpcb_vec_init(3 * len);
    u = t + len;
    v = u + len;

    /* sin, cos of h0 */
    fmpcb_sin_cos(s0, c0, h, prec);

    /* t = tan((h-h0)/2) */
    fmpcb_zero(u);
    _fmpcb_vec_scalar_mul_2exp_si(u + 1, h + 1, hlen - 1, -1);
    _fmpcb_poly_tan_series(t, u, hlen, len, prec);

    /* v = 1 + t^2 */
    _fmpcb_poly_mullow(v, t, len, t, len, len, prec);
    fmpcb_add_ui(v, v, 1, prec);

    /* u = 1/(1+t^2) */
    _fmpcb_poly_inv_series(u, v, len, len, prec);

    /* sine */
    _fmpcb_poly_mullow(s, t, len, u, len, len, prec);
    _fmpcb_vec_scalar_mul_2exp_si(s, s, len, 1);

    /* cosine */
    fmpcb_sub_ui(v, v, 2, prec);
    _fmpcb_vec_neg(v, v, len);
    _fmpcb_poly_mullow(c, v, len, u, len, len, prec);

    /* sin(h0 + h1) = cos(h0) sin(h1) + sin(h0) cos(h1)
       cos(h0 + h1) = cos(h0) cos(h1) - sin(h0) sin(h1) */
    if (!fmpcb_is_zero(s0))
    {
        _fmpcb_vec_scalar_mul(t, s, len, c0, prec);
        _fmpcb_vec_scalar_mul(u, c, len, s0, prec);
        _fmpcb_vec_scalar_mul(v, s, len, s0, prec);
        _fmpcb_vec_add(s, t, u, len, prec);
        _fmpcb_vec_scalar_mul(t, c, len, c0, prec);
        _fmpcb_vec_sub(c, t, v, len, prec);
    }

    _fmpcb_vec_clear(t, 3 * len);

    fmpcb_clear(s0);
    fmpcb_clear(c0);
}

void
fmpcb_poly_sin_cos_series_tangent(fmpcb_poly_t s, fmpcb_poly_t c,
                                    const fmpcb_poly_t h, long n, long prec)
{
    long hlen = h->length;

    if (n == 0)
    {
        fmpcb_poly_zero(s);
        fmpcb_poly_zero(c);
        return;
    }

    if (hlen == 0)
    {
        fmpcb_poly_zero(s);
        fmpcb_poly_one(c);
        return;
    }

    fmpcb_poly_fit_length(s, n);
    fmpcb_poly_fit_length(c, n);
    _fmpcb_poly_sin_cos_series_tangent(s->coeffs, c->coeffs,
        h->coeffs, hlen, n, prec);
    _fmpcb_poly_set_length(s, n);
    _fmpcb_poly_normalise(s);
    _fmpcb_poly_set_length(c, n);
    _fmpcb_poly_normalise(c);
}

