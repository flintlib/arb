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

#define MULLOW(z, x, xn, y, yn, nn, prec) \
    if ((xn) >= (yn)) \
        _fmprb_poly_mullow(z, x, xn, y, yn, nn, prec); \
    else \
        _fmprb_poly_mullow(z, y, yn, x, xn, nn, prec); \

void
_fmprb_poly_inv_series(fmprb_ptr Qinv,
    fmprb_srcptr Q, long Qlen, long len, long prec)
{
    fmprb_ui_div(Qinv, 1UL, Q, prec);

    if (Qlen == 1)
    {
        _fmprb_vec_zero(Qinv + 1, len - 1);
    }
    else
    {
        long Qnlen, Wlen, W2len;
        fmprb_ptr W;

        W = _fmprb_vec_init(len);

        NEWTON_INIT(1, len)
        NEWTON_LOOP(m, n)

        Qnlen = FLINT_MIN(Qlen, n);
        Wlen = FLINT_MIN(Qnlen + m - 1, n);
        W2len = Wlen - m;
        MULLOW(W, Q, Qnlen, Qinv, m, Wlen, prec);
        MULLOW(Qinv + m, Qinv, m, W + m, W2len, n - m, prec);
        _fmprb_vec_neg(Qinv + m, Qinv + m, n - m);

        NEWTON_END_LOOP
        NEWTON_END

        _fmprb_vec_clear(W, len);
    }
}

void
fmprb_poly_inv_series(fmprb_poly_t Qinv, const fmprb_poly_t Q, long n, long prec)
{
    if (n == 0 || Q->length == 0)
    {
        printf("fmprb_poly_inv_series: require n > 0 and nonzero input\n");
        abort();
    }

    if (Qinv == Q)
    {
        fmprb_poly_t t;
        fmprb_poly_init(t);
        fmprb_poly_inv_series(t, Q, n, prec);
        fmprb_poly_swap(Qinv, t);
        fmprb_poly_clear(t);
        return;
    }

    fmprb_poly_fit_length(Qinv, n);
    _fmprb_poly_inv_series(Qinv->coeffs, Q->coeffs, Q->length, n, prec);
    _fmprb_poly_set_length(Qinv, n);
    _fmprb_poly_normalise(Qinv);
}

