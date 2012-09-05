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
_fmprb_poly_neg(fmprb_struct * B, const fmprb_struct * A, long n, long prec)
{
    long i;

    for (i = 0; i < n; i++)
        fmprb_neg(B + i, A + i);
}

void
_fmprb_poly_inv_series(fmprb_struct * Qinv, const fmprb_struct * Q, long len, long prec)
{
    long a[FLINT_BITS];
    long i, m, n;
    fmprb_struct * W;

    W = flint_malloc(len * sizeof(fmprb_struct));
    for (i = 0; i < len; i++)
        fmprb_init(W + i);

    a[i = 0] = n = len;
    while (n >= 2)
        a[++i] = (n = (n + 1) / 2);

    fmprb_ui_div(Qinv, 1UL, Q, prec);

    for (i--; i >= 0; i--)
    {
        m = n;
        n = a[i];

        _fmprb_poly_mullow(W, Q, n, Qinv, m, n, prec);
        _fmprb_poly_mullow(Qinv + m, Qinv, m, W + m, n - m, n - m, prec);
        _fmprb_poly_neg(Qinv + m, Qinv + m, n - m, prec);
    }

    for (i = 0; i < len; i++)
        fmprb_clear(W + i);

    flint_free(W);
}

void
fmprb_poly_inv_series(fmprb_poly_t Qinv, const fmprb_poly_t Q, long n, long prec)
{
    if (n == 0 || Q->length == 0)
        abort();

    if (Q->length < n || Qinv == Q)
    {
        long i;
        fmprb_poly_t t;

        fmprb_poly_init(t);
        fmprb_poly_fit_length(t, n);

        for (i = 0; i < FLINT_MIN(n, Q->length); i++)
            fmprb_set(t->coeffs + i, Q->coeffs + i);

        for (i = FLINT_MIN(n, Q->length); i < n; i++)
            fmprb_zero(t->coeffs + i);

        _fmprb_poly_set_length(t, n);
        fmprb_poly_inv_series(Qinv, t, n, prec);

        fmprb_poly_clear(t);
        return;
    }

    fmprb_poly_fit_length(Qinv, n);
    _fmprb_poly_inv_series(Qinv->coeffs, Q->coeffs, n, prec);
    _fmprb_poly_set_length(Qinv, n);
    _fmprb_poly_normalise(Qinv);
}
