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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "fmpcb_poly.h"

void
_fmpcb_poly_exp_series_basecase(fmpcb_ptr f,
        fmpcb_srcptr h, long hlen, long n, long prec)
{
    long j, k, alen = FLINT_MIN(n, hlen);
    fmpcb_ptr a;
    fmpcb_t s;

    fmpcb_init(s);
    a = _fmpcb_vec_init(alen);

    fmpcb_exp(f, h, prec);

    for (k = 1; k < alen; k++)
        fmpcb_mul_ui(a + k, h + k, k, prec);

    for (k = 1; k < n; k++)
    {
        fmpcb_zero(s);
        for (j = 1; j < FLINT_MIN(k + 1, hlen); j++)
            fmpcb_addmul(s, a + j, f + k - j, prec);

        fmpcb_div_ui(f + k, s, k, prec);
    }

    fmpcb_clear(s);
    _fmpcb_vec_clear(a, alen);
}

void
fmpcb_poly_exp_series_basecase(fmpcb_poly_t f, const fmpcb_poly_t h, long n, long prec)
{
    long hlen = h->length;

    if (n == 0)
    {
        fmpcb_poly_zero(f);
        return;
    }

    if (hlen == 0)
    {
        fmpcb_poly_one(f);
        return;
    }

    fmpcb_poly_fit_length(f, n);
    _fmpcb_poly_exp_series_basecase(f->coeffs, h->coeffs, hlen, n, prec);
    _fmpcb_poly_set_length(f, n);
    _fmpcb_poly_normalise(f);
}
