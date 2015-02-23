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

#include "acb_poly.h"

#define TANGENT_CUTOFF 80

void
_acb_poly_sin_cos_series(acb_ptr s, acb_ptr c, const acb_srcptr h, long hlen, long n, long prec)
{
    hlen = FLINT_MIN(hlen, n);

    if (hlen == 1)
    {
        acb_sin_cos(s, c, h, prec);
        _acb_vec_zero(s + 1, n - 1);
        _acb_vec_zero(c + 1, n - 1);
    }
    else if (n == 2)
    {
        acb_t t;
        acb_init(t);
        acb_set(t, h + 1);
        acb_sin_cos(s, c, h, prec);
        acb_mul(s + 1, c, t, prec);
        acb_neg(t, t);
        acb_mul(c + 1, s, t, prec);
        acb_clear(t);
    }
    else if (hlen < TANGENT_CUTOFF)
        _acb_poly_sin_cos_series_basecase(s, c, h, hlen, n, prec, 0);
    else
        _acb_poly_sin_cos_series_tangent(s, c, h, hlen, n, prec, 0);
}

void
acb_poly_sin_cos_series(acb_poly_t s, acb_poly_t c,
                                    const acb_poly_t h, long n, long prec)
{
    long hlen = h->length;

    if (n == 0)
    {
        acb_poly_zero(s);
        acb_poly_zero(c);
        return;
    }

    if (hlen == 0)
    {
        acb_poly_zero(s);
        acb_poly_one(c);
        return;
    }

    if (hlen == 1)
        n = 1;

    acb_poly_fit_length(s, n);
    acb_poly_fit_length(c, n);
    _acb_poly_sin_cos_series(s->coeffs, c->coeffs, h->coeffs, hlen, n, prec);
    _acb_poly_set_length(s, n);
    _acb_poly_normalise(s);
    _acb_poly_set_length(c, n);
    _acb_poly_normalise(c);
}

