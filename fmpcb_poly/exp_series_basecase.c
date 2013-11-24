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

#define MUL_CUTOFF 24

static void
_fmpcb_poly_exp_series_basecase_rec(fmpcb_ptr f, fmpcb_ptr a,
        fmpcb_srcptr h, long hlen, long n, long prec)
{
    long j, k;

    fmpcb_t s;
    fmpcb_init(s);

    fmpcb_exp(f, h, prec);

    for (k = 1; k < hlen; k++)
        fmpcb_mul_ui(a + k, h + k, k, prec);

    for (k = 1; k < n; k++)
    {
        fmpcb_zero(s);
        for (j = 1; j < FLINT_MIN(k + 1, hlen); j++)
            fmpcb_addmul(s, a + j, f + k - j, prec);

        fmpcb_div_ui(f + k, s, k, prec);
    }

    fmpcb_clear(s);
}

void
_fmpcb_poly_exp_series_basecase(fmpcb_ptr f,
        fmpcb_srcptr h, long hlen, long n, long prec)
{
    hlen = FLINT_MIN(n, hlen);

    if (n < MUL_CUTOFF || hlen < 0.9 * n)
    {
        fmpcb_ptr t = _fmpcb_vec_init(hlen);
        _fmpcb_poly_exp_series_basecase_rec(f, t, h, hlen, n, prec);
        _fmpcb_vec_clear(t, hlen);
    }
    else
    {
        long m, v;
        fmpcb_ptr t, u;

        m = (n + 2) / 3;
        v = m * 2;

        t = _fmpcb_vec_init(n);
        u = _fmpcb_vec_init(n - m);

        _fmpcb_poly_mullow(t, h + m, hlen - m, h + m, hlen - m, n - v, prec);
        _fmpcb_vec_scalar_mul_2exp_si(t, t, n - v, -1);

        _fmpcb_vec_set(u, h + m, v - m);
        _fmpcb_poly_add(u + v - m, t, n - v, h + v, hlen - v, prec);
        _fmpcb_poly_exp_series_basecase_rec(f, t, h, m, n, prec);
        _fmpcb_poly_mullow(t, f, n, u, n - m, n - m, prec);
        _fmpcb_poly_add(f + m, f + m, n - m, t, n - m, prec);

        _fmpcb_vec_clear(t, n);
        _fmpcb_vec_clear(u, n - m);
    }
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
