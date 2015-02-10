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

void
_acb_poly_sin_cos_series_basecase(acb_ptr s, acb_ptr c, acb_srcptr h, long hlen,
        long n, long prec, int times_pi)
{
    long j, k, alen = FLINT_MIN(n, hlen);
    acb_ptr a;
    acb_t t, u;

    if (times_pi)
        acb_sin_cos_pi(s, c, h, prec);
    else
        acb_sin_cos(s, c, h, prec);

    if (hlen == 1)
    {
        _acb_vec_zero(s + 1, n - 1);
        _acb_vec_zero(c + 1, n - 1);
        return;
    }

    acb_init(t);
    acb_init(u);
    a = _acb_vec_init(alen);

    for (k = 1; k < alen; k++)
        acb_mul_ui(a + k, h + k, k, prec);

    if (times_pi)
    {
        acb_const_pi(t, prec);
        _acb_vec_scalar_mul(a + 1, a + 1, alen - 1, t, prec);
    }

    for (k = 1; k < n; k++)
    {
        acb_zero(t);
        acb_zero(u);

        for (j = 1; j < FLINT_MIN(k + 1, hlen); j++)
        {
            acb_submul(t, a + j, s + k - j, prec);
            acb_addmul(u, a + j, c + k - j, prec);
        }

        acb_div_ui(c + k, t, k, prec);
        acb_div_ui(s + k, u, k, prec);
    }

    acb_clear(t);
    acb_clear(u);
    _acb_vec_clear(a, alen);
}

void
acb_poly_sin_cos_series_basecase(acb_poly_t s, acb_poly_t c,
        const acb_poly_t h, long n, long prec, int times_pi)
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

    acb_poly_fit_length(s, n);
    acb_poly_fit_length(c, n);
    _acb_poly_sin_cos_series_basecase(s->coeffs, c->coeffs, h->coeffs, hlen, n, prec, times_pi);
    _acb_poly_set_length(s, n);
    _acb_poly_normalise(s);
    _acb_poly_set_length(c, n);
    _acb_poly_normalise(c);
}

