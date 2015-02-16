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

#include "acb_poly.h"

#define MULLOW(z, x, xn, y, yn, nn, prec) \
    if ((xn) >= (yn)) \
        _acb_poly_mullow(z, x, xn, y, yn, nn, prec); \
    else \
        _acb_poly_mullow(z, y, yn, x, xn, nn, prec); \

void
_acb_poly_inv_series(acb_ptr Qinv,
    acb_srcptr Q, long Qlen, long len, long prec)
{
    acb_inv(Qinv, Q, prec);

    if (Qlen == 1)
    {
        _acb_vec_zero(Qinv + 1, len - 1);
    }
    else if (len == 2)
    {
        acb_div(Qinv + 1, Qinv, Q, prec);
        acb_mul(Qinv + 1, Qinv + 1, Q + 1, prec);
        acb_neg(Qinv + 1, Qinv + 1);
    }
    else
    {
        long Qnlen, Wlen, W2len;
        acb_ptr W;

        W = _acb_vec_init(len);

        NEWTON_INIT(1, len)
        NEWTON_LOOP(m, n)

        Qnlen = FLINT_MIN(Qlen, n);
        Wlen = FLINT_MIN(Qnlen + m - 1, n);
        W2len = Wlen - m;
        MULLOW(W, Q, Qnlen, Qinv, m, Wlen, prec);
        MULLOW(Qinv + m, Qinv, m, W + m, W2len, n - m, prec);
        _acb_vec_neg(Qinv + m, Qinv + m, n - m);

        NEWTON_END_LOOP
        NEWTON_END

        _acb_vec_clear(W, len);
    }
}

void
acb_poly_inv_series(acb_poly_t Qinv, const acb_poly_t Q, long n, long prec)
{
    if (n == 0)
    {
        acb_poly_zero(Qinv);
        return;
    }

    if (Q->length == 0)
    {
        acb_poly_fit_length(Qinv, n);
        _acb_vec_indeterminate(Qinv->coeffs, n);
        _acb_poly_set_length(Qinv, n);
        return;
    }

    if (Qinv == Q)
    {
        acb_poly_t t;
        acb_poly_init(t);
        acb_poly_inv_series(t, Q, n, prec);
        acb_poly_swap(Qinv, t);
        acb_poly_clear(t);
        return;
    }

    acb_poly_fit_length(Qinv, n);
    _acb_poly_inv_series(Qinv->coeffs, Q->coeffs, Q->length, n, prec);
    _acb_poly_set_length(Qinv, n);
    _acb_poly_normalise(Qinv);
}

