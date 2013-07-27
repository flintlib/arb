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
_fmpcb_poly_sin_cos_series_basecase(fmpcb_ptr s,
                                    fmpcb_ptr c, fmpcb_srcptr h, long hlen, long n, long prec)
{
    long j, k, alen = FLINT_MIN(n, hlen);
    fmpcb_ptr a;
    fmpcb_t t, u;

    fmpcb_sin_cos(s, c, h, prec);

    if (hlen == 1)
    {
        _fmpcb_vec_zero(s + 1, n - 1);
        _fmpcb_vec_zero(c + 1, n - 1);
        return;
    }

    fmpcb_init(t);
    fmpcb_init(u);
    a = _fmpcb_vec_init(alen);

    for (k = 1; k < alen; k++)
        fmpcb_mul_ui(a + k, h + k, k, prec);

    for (k = 1; k < n; k++)
    {
        fmpcb_zero(t);
        fmpcb_zero(u);

        for (j = 1; j < FLINT_MIN(k + 1, hlen); j++)
        {
            fmpcb_submul(t, a + j, s + k - j, prec);
            fmpcb_addmul(u, a + j, c + k - j, prec);
        }

        fmpcb_div_ui(c + k, t, k, prec);
        fmpcb_div_ui(s + k, u, k, prec);
    }

    fmpcb_clear(t);
    fmpcb_clear(u);
    _fmpcb_vec_clear(a, alen);
}

void
fmpcb_poly_sin_cos_series_basecase(fmpcb_poly_t s, fmpcb_poly_t c,
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
    _fmpcb_poly_sin_cos_series_basecase(s->coeffs, c->coeffs, h->coeffs, hlen, n, prec);
    _fmpcb_poly_set_length(s, n);
    _fmpcb_poly_normalise(s);
    _fmpcb_poly_set_length(c, n);
    _fmpcb_poly_normalise(c);
}

