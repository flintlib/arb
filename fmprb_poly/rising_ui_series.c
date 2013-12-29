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

    Copyright (C) 2012, 2013 Fredrik Johansson

******************************************************************************/

#include "fmprb_poly.h"
#include "gamma.h"

static void
_fmprb_poly_rising_ui_series_bsplit(fmprb_ptr res,
    fmprb_srcptr f, long flen, ulong a, ulong b,
        long trunc, long prec)
{
    flen = FLINT_MIN(flen, trunc);

    if (b - a == 1)
    {
        fmprb_add_ui(res, f, a, prec);
        _fmprb_vec_set(res + 1, f + 1, flen - 1);
    }
    else
    {
        fmprb_ptr L, R;
        long len1, len2;

        long m = a + (b - a) / 2;

        len1 = poly_pow_length(flen, m - a, trunc);
        len2 = poly_pow_length(flen, b - m, trunc);

        L = _fmprb_vec_init(len1 + len2);
        R = L + len1;

        _fmprb_poly_rising_ui_series_bsplit(L, f, flen, a, m, trunc, prec);
        _fmprb_poly_rising_ui_series_bsplit(R, f, flen, m, b, trunc, prec);

        _fmprb_poly_mullow(res, L, len1, R, len2,
            FLINT_MIN(trunc, len1 + len2 - 1), prec);

        _fmprb_vec_clear(L, len1 + len2);
    }
}

void
_fmprb_poly_rising_ui_series(fmprb_ptr res,
    fmprb_srcptr f, long flen, ulong r,
        long trunc, long prec)
{
    if (trunc == 1 || flen == 1)
    {
        gamma_rising_fmprb_ui_bsplit(res, f, r, prec);
        _fmprb_vec_zero(res + 1, trunc - 1);
    }
    else if (trunc == 2)
    {
        gamma_rising2_fmprb_ui(res, res + 1, f, r, prec);
        fmprb_mul(res + 1, res + 1, f + 1, prec);
    }
    else
    {
        _fmprb_poly_rising_ui_series_bsplit(res, f, flen, 0, r, trunc, prec);
    }
}

void
fmprb_poly_rising_ui_series(fmprb_poly_t res, const fmprb_poly_t f, ulong r, long trunc, long prec)
{
    long len;

    if (f->length == 0 && r != 0)
    {
        fmprb_poly_zero(res);
        return;
    }

    if (r == 0)
    {
        fmprb_poly_one(res);
        return;
    }

    len = poly_pow_length(f->length, r, trunc);

    if (f == res)
    {
        fmprb_poly_t tmp;
        fmprb_poly_init(tmp);
        fmprb_poly_rising_ui_series(tmp, f, r, len, prec);
        fmprb_poly_swap(tmp, res);
        fmprb_poly_clear(tmp);
    }
    else
    {
        fmprb_poly_fit_length(res, len);
        _fmprb_poly_rising_ui_series(res->coeffs, f->coeffs, f->length, r, len, prec);
        _fmprb_poly_set_length(res, len);
        _fmprb_poly_normalise(res);
    }
}

