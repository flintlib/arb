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

#include "fmprb_poly.h"

#define MUL(z, zlen, x, xlen, y, ylen, trunc, prec) \
    do { \
        long slen = FLINT_MIN(xlen + ylen - 1, trunc); \
        _fmprb_poly_mullow(z, x, xlen, y, ylen, slen, prec); \
        zlen = slen; \
    } while (0)

void
_fmprb_poly_pow_ui_trunc_binexp(fmprb_ptr res,
    fmprb_srcptr f, long flen, ulong exp, long len, long prec)
{
    fmprb_ptr v, R, S, T;
    long rlen;
    ulong bit;

    if (exp <= 1)
    {
        if (exp == 0)
            fmprb_one(res);
        else if (exp == 1)
            _fmprb_vec_set_round(res, f, len, prec);
        return;
    }

    /* (f * x^r)^m = x^(rm) * f^m */
    while (flen > 1 && fmprb_is_zero(f))
    {
        if (((ulong) len) > exp)
        {
            _fmprb_vec_zero(res, exp);
            len -= exp;
            res += exp;
        }
        else
        {
            _fmprb_vec_zero(res, len);
            return;
        }

        f++;
        flen--;
    }

    if (exp == 2)
    {
        _fmprb_poly_mullow(res, f, flen, f, flen, len, prec);
        return;
    }

    if (flen == 1)
    {
        fmprb_pow_ui(res, f, exp, prec);
        return;
    }

    v = _fmprb_vec_init(len);
    bit = 1UL << (FLINT_BIT_COUNT(exp) - 2);
    
    if (n_zerobits(exp) % 2)
    {
        R = res;
        S = v;
    }
    else
    {
        R = v;
        S = res;
    }

    MUL(R, rlen, f, flen, f, flen, len, prec);

    if (bit & exp)
    {
        MUL(S, rlen, R, rlen, f, flen, len, prec);
        T = R;
        R = S;
        S = T;
    }
    
    while (bit >>= 1)
    {
        if (bit & exp)
        {
            MUL(S, rlen, R, rlen, R, rlen, len, prec);
            MUL(R, rlen, S, rlen, f, flen, len, prec);
        }
        else
        {
            MUL(S, rlen, R, rlen, R, rlen, len, prec);
            T = R;
            R = S;
            S = T;
        }
    }
    
    _fmprb_vec_clear(v, len);
}

void
fmprb_poly_pow_ui_trunc_binexp(fmprb_poly_t res,
    const fmprb_poly_t poly, ulong exp, long len, long prec)
{
    long flen, rlen;

    flen = poly->length;

    if (exp == 0 && len != 0)
    {
        fmprb_poly_one(res);
    }
    else if (flen == 0 || len == 0)
    {
        fmprb_poly_zero(res);
    }
    else
    {
        rlen = poly_pow_length(flen, exp, len);

        if (res != poly)
        {
            fmprb_poly_fit_length(res, rlen);
            _fmprb_poly_pow_ui_trunc_binexp(res->coeffs,
                poly->coeffs, flen, exp, rlen, prec);
            _fmprb_poly_set_length(res, rlen);
            _fmprb_poly_normalise(res);
        }
        else
        {
            fmprb_poly_t t;
            fmprb_poly_init2(t, rlen);
            _fmprb_poly_pow_ui_trunc_binexp(t->coeffs,
                poly->coeffs, flen, exp, rlen, prec);
            _fmprb_poly_set_length(t, rlen);
            _fmprb_poly_normalise(t);
            fmprb_poly_swap(res, t);
            fmprb_poly_clear(t);
        }
    }
}

