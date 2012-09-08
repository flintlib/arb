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

#include "fmprb_poly.h"

void
_fmprb_poly_exp_series_basecase(fmprb_struct * f,
        const fmprb_struct * h, long hlen, long n, long prec)
{
    long j, k, alen = FLINT_MIN(n, hlen);
    fmprb_struct * a;
    fmprb_t s;

    fmprb_init(s);
    a = _fmprb_vec_init(alen);

    fmprb_exp(f, h, prec);

    for (k = 1; k < alen; k++)
        fmprb_mul_ui(a + k, h + k, k, prec);

    for (k = 1; k < n; k++)
    {
        fmprb_mul(s, a + 1, f + k - 1, prec);
        for (j = 2; j < FLINT_MIN(k + 1, hlen); j++)
            fmprb_addmul(s, a + j, f + k - j, prec);
        fmprb_div_ui(f + k, s, k, prec);
    }

    fmprb_clear(s);
    _fmprb_vec_clear(a, alen);
}

void
fmprb_poly_exp_series_basecase(fmprb_poly_t f, const fmprb_poly_t h, long n, long prec)
{
    long hlen = h->length;

    if (n == 0)
    {
        fmprb_poly_zero(f);
        return;
    }

    if (hlen == 0)
    {
        fmprb_poly_one(f);
        return;
    }

    fmprb_poly_fit_length(f, n);
    _fmprb_poly_exp_series_basecase(f->coeffs, h->coeffs, hlen, n, prec);
    _fmprb_poly_set_length(f, n);
    _fmprb_poly_normalise(f);
}
