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

#define MUL_CUTOFF 24

static void
_fmprb_poly_exp_series_basecase_rec(fmprb_ptr f, fmprb_ptr a,
        fmprb_srcptr h, long hlen, long n, long prec)
{
    long j, k;

    fmprb_t s;
    fmprb_init(s);

    fmprb_exp(f, h, prec);

    for (k = 1; k < hlen; k++)
        fmprb_mul_ui(a + k, h + k, k, prec);

    for (k = 1; k < n; k++)
    {
        fmprb_zero(s);
        for (j = 1; j < FLINT_MIN(k + 1, hlen); j++)
            fmprb_addmul(s, a + j, f + k - j, prec);

        fmprb_div_ui(f + k, s, k, prec);
    }

    fmprb_clear(s);
}

void
_fmprb_poly_exp_series_basecase(fmprb_ptr f,
        fmprb_srcptr h, long hlen, long n, long prec)
{
    hlen = FLINT_MIN(n, hlen);

    if (n < MUL_CUTOFF || hlen < 0.9 * n)
    {
        fmprb_ptr t = _fmprb_vec_init(hlen);
        _fmprb_poly_exp_series_basecase_rec(f, t, h, hlen, n, prec);
        _fmprb_vec_clear(t, hlen);
    }
    else
    {
        long m, v;
        fmprb_ptr t, u;

        m = (n + 2) / 3;
        v = m * 2;

        t = _fmprb_vec_init(n);
        u = _fmprb_vec_init(n - m);

        _fmprb_poly_mullow(t, h + m, hlen - m, h + m, hlen - m, n - v, prec);
        _fmprb_vec_scalar_mul_2exp_si(t, t, n - v, -1);

        _fmprb_vec_set(u, h + m, v - m);
        _fmprb_poly_add(u + v - m, t, n - v, h + v, hlen - v, prec);
        _fmprb_poly_exp_series_basecase_rec(f, t, h, m, n, prec);
        _fmprb_poly_mullow(t, f, n, u, n - m, n - m, prec);
        _fmprb_poly_add(f + m, f + m, n - m, t, n - m, prec);

        _fmprb_vec_clear(t, n);
        _fmprb_vec_clear(u, n - m);
    }
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
