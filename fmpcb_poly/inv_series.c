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
_fmpcb_poly_inv_series(fmpcb_struct * Qinv, const fmpcb_struct * Q, long len, long prec)
{
    long a[FLINT_BITS];
    long i, m, n;
    fmpcb_struct * W;

    W = _fmpcb_vec_init(len);

    a[i = 0] = n = len;
    while (n >= 2)
        a[++i] = (n = (n + 1) / 2);

    fmpcb_inv(Qinv, Q, prec);

    for (i--; i >= 0; i--)
    {
        m = n;
        n = a[i];

        _fmpcb_poly_mullow(W, Q, n, Qinv, m, n, prec);
        _fmpcb_poly_mullow(Qinv + m, Qinv, m, W + m, n - m, n - m, prec);
        _fmpcb_vec_neg(Qinv + m, Qinv + m, n - m);
    }

    _fmpcb_vec_clear(W, len);
}

void
fmpcb_poly_inv_series(fmpcb_poly_t Qinv, const fmpcb_poly_t Q, long n, long prec)
{
    if (n == 0 || Q->length == 0)
    {
        printf("fmpcb_poly_inv_series: input must be nonzero\n");
        abort();
    }

    if (Q->length < n || Qinv == Q)
    {
        long i;
        fmpcb_poly_t t;

        fmpcb_poly_init(t);
        fmpcb_poly_fit_length(t, n);

        for (i = 0; i < FLINT_MIN(n, Q->length); i++)
            fmpcb_set(t->coeffs + i, Q->coeffs + i);

        for (i = FLINT_MIN(n, Q->length); i < n; i++)
            fmpcb_zero(t->coeffs + i);

        _fmpcb_poly_set_length(t, n);
        fmpcb_poly_inv_series(Qinv, t, n, prec);

        fmpcb_poly_clear(t);
        return;
    }

    fmpcb_poly_fit_length(Qinv, n);
    _fmpcb_poly_inv_series(Qinv->coeffs, Q->coeffs, n, prec);
    _fmpcb_poly_set_length(Qinv, n);
    _fmpcb_poly_normalise(Qinv);
}

