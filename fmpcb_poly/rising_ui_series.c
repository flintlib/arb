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

#include "fmpcb_poly.h"
#include "gamma.h"

static void
_fmpcb_poly_rising_ui_series_bsplit(fmpcb_ptr res,
    fmpcb_srcptr f, long flen, ulong a, ulong b,
        long trunc, long prec)
{
    flen = FLINT_MIN(flen, trunc);

    if (b - a == 1)
    {
        fmpcb_add_ui(res, f, a, prec);
        _fmpcb_vec_set(res + 1, f + 1, flen - 1);
    }
    else
    {
        fmpcb_ptr L, R;
        long len1, len2;

        long m = a + (b - a) / 2;

        len1 = poly_pow_length(flen, m - a, trunc);
        len2 = poly_pow_length(flen, b - m, trunc);

        L = _fmpcb_vec_init(len1 + len2);
        R = L + len1;

        _fmpcb_poly_rising_ui_series_bsplit(L, f, flen, a, m, trunc, prec);
        _fmpcb_poly_rising_ui_series_bsplit(R, f, flen, m, b, trunc, prec);

        _fmpcb_poly_mullow(res, L, len1, R, len2,
            FLINT_MIN(trunc, len1 + len2 - 1), prec);

        _fmpcb_vec_clear(L, len1 + len2);
    }
}

void
_fmpcb_poly_rising_ui_series(fmpcb_ptr res,
    fmpcb_srcptr f, long flen, ulong r,
        long trunc, long prec)
{
    if (trunc == 1 || flen == 1)
    {
        gamma_rising_fmpcb_ui_bsplit(res, f, r, prec);
        _fmpcb_vec_zero(res + 1, trunc - 1);
    }
    else if (trunc == 2)
    {
        gamma_rising2_fmpcb_ui(res, res + 1, f, r, prec);
        fmpcb_mul(res + 1, res + 1, f + 1, prec);
    }
    else
    {
        _fmpcb_poly_rising_ui_series_bsplit(res, f, flen, 0, r, trunc, prec);
    }
}

void
fmpcb_poly_rising_ui_series(fmpcb_poly_t res, const fmpcb_poly_t f, ulong r, long trunc, long prec)
{
    long len;

    if (f->length == 0 && r != 0)
    {
        fmpcb_poly_zero(res);
        return;
    }

    if (r == 0)
    {
        fmpcb_poly_one(res);
        return;
    }

    len = poly_pow_length(f->length, r, trunc);

    if (f == res)
    {
        fmpcb_poly_t tmp;
        fmpcb_poly_init(tmp);
        fmpcb_poly_rising_ui_series(tmp, f, r, len, prec);
        fmpcb_poly_swap(tmp, res);
        fmpcb_poly_clear(tmp);
    }
    else
    {
        fmpcb_poly_fit_length(res, len);
        _fmpcb_poly_rising_ui_series(res->coeffs, f->coeffs, f->length, r, len, prec);
        _fmpcb_poly_set_length(res, len);
        _fmpcb_poly_normalise(res);
    }
}

