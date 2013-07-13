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
_fmprb_poly_div_series(fmprb_struct * Q, const fmprb_struct * A, long Alen,
    const fmprb_struct * B, long Blen, long n, long prec)
{
    Alen = FLINT_MIN(Alen, n);
    Blen = FLINT_MIN(Blen, n);

    if (Blen == 1)
    {
        _fmprb_vec_scalar_div(Q, A, Alen, B, prec);
        _fmprb_vec_zero(Q + Alen, n - Alen);
    }
    else
    {
        fmprb_struct * Binv;
        Binv = _fmprb_vec_init(n);
        _fmprb_poly_inv_series(Binv, B, Blen, n, prec);
        MULLOW(Q, A, Alen, Binv, n, n, prec);
        _fmprb_vec_clear(Binv, n);
    }
}

void
fmprb_poly_div_series(fmprb_poly_t Q, const fmprb_poly_t A, const fmprb_poly_t B, long n, long prec)
{
    if (n == 0 || B->length == 0)
    {
        printf("fmprb_poly_inv_series: require n > 0 and nonzero input\n");
        abort();
    }

    if (A->length == 0)
    {
        fmprb_poly_zero(Q);
        return;
    }

    if (Q == A || Q == B)
    {
        fmprb_poly_t t;
        fmprb_poly_init(t);
        fmprb_poly_div_series(t, A, B, n, prec);
        fmprb_poly_swap(Q, t);
        fmprb_poly_clear(t);
        return;
    }

    fmprb_poly_fit_length(Q, n);
    _fmprb_poly_div_series(Q->coeffs, A->coeffs, A->length, B->coeffs, B->length, n, prec);
    _fmprb_poly_set_length(Q, n);
    _fmprb_poly_normalise(Q);
}

