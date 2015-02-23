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

#include "arb_poly.h"

static void
_arb_poly_rising_ui_series_bsplit(arb_ptr res,
    arb_srcptr f, long flen, ulong a, ulong b,
        long trunc, long prec)
{
    flen = FLINT_MIN(flen, trunc);

    if (b - a == 1)
    {
        arb_add_ui(res, f, a, prec);
        _arb_vec_set(res + 1, f + 1, flen - 1);
    }
    else
    {
        arb_ptr L, R;
        long len1, len2;

        long m = a + (b - a) / 2;

        len1 = poly_pow_length(flen, m - a, trunc);
        len2 = poly_pow_length(flen, b - m, trunc);

        L = _arb_vec_init(len1 + len2);
        R = L + len1;

        _arb_poly_rising_ui_series_bsplit(L, f, flen, a, m, trunc, prec);
        _arb_poly_rising_ui_series_bsplit(R, f, flen, m, b, trunc, prec);

        _arb_poly_mullow(res, L, len1, R, len2,
            FLINT_MIN(trunc, len1 + len2 - 1), prec);

        _arb_vec_clear(L, len1 + len2);
    }
}

void
_arb_poly_rising_ui_series(arb_ptr res,
    arb_srcptr f, long flen, ulong r,
        long trunc, long prec)
{
    if (trunc == 1 || flen == 1)
    {
        arb_rising_ui(res, f, r, prec);
        _arb_vec_zero(res + 1, trunc - 1);
    }
    else if (trunc == 2)
    {
        arb_rising2_ui(res, res + 1, f, r, prec);
        arb_mul(res + 1, res + 1, f + 1, prec);
    }
    else
    {
        _arb_poly_rising_ui_series_bsplit(res, f, flen, 0, r, trunc, prec);
    }
}

void
arb_poly_rising_ui_series(arb_poly_t res, const arb_poly_t f, ulong r, long trunc, long prec)
{
    long len;

    if ((f->length == 0 && r != 0) || trunc == 0)
    {
        arb_poly_zero(res);
        return;
    }

    if (r == 0)
    {
        arb_poly_one(res);
        return;
    }

    len = poly_pow_length(f->length, r, trunc);

    if (f == res)
    {
        arb_poly_t tmp;
        arb_poly_init(tmp);
        arb_poly_rising_ui_series(tmp, f, r, len, prec);
        arb_poly_swap(tmp, res);
        arb_poly_clear(tmp);
    }
    else
    {
        arb_poly_fit_length(res, len);
        _arb_poly_rising_ui_series(res->coeffs, f->coeffs, f->length, r, len, prec);
        _arb_poly_set_length(res, len);
        _arb_poly_normalise(res);
    }
}

