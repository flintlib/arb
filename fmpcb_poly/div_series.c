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

void 
_fmpcb_poly_div_series(fmpcb_ptr Q, fmpcb_srcptr A, long Alen,
    fmpcb_srcptr B, long Blen, long n, long prec)
{
    Alen = FLINT_MIN(Alen, n);
    Blen = FLINT_MIN(Blen, n);

    if (Blen == 1)
    {
        _fmpcb_vec_scalar_div(Q, A, Alen, B, prec);
        _fmpcb_vec_zero(Q + Alen, n - Alen);
    }
    else
    {
        fmpcb_ptr Binv;
        Binv = _fmpcb_vec_init(n);
        _fmpcb_poly_inv_series(Binv, B, Blen, n, prec);
        _fmpcb_poly_mullow(Q, Binv, n, A, Alen, n, prec);
        _fmpcb_vec_clear(Binv, n);
    }
}

void
fmpcb_poly_div_series(fmpcb_poly_t Q, const fmpcb_poly_t A, const fmpcb_poly_t B, long n, long prec)
{
    if (n == 0 || B->length == 0)
    {
        printf("fmpcb_poly_inv_series: require n > 0 and nonzero input\n");
        abort();
    }

    if (A->length == 0)
    {
        fmpcb_poly_zero(Q);
        return;
    }

    if (Q == A || Q == B)
    {
        fmpcb_poly_t t;
        fmpcb_poly_init(t);
        fmpcb_poly_div_series(t, A, B, n, prec);
        fmpcb_poly_swap(Q, t);
        fmpcb_poly_clear(t);
        return;
    }

    fmpcb_poly_fit_length(Q, n);
    _fmpcb_poly_div_series(Q->coeffs, A->coeffs, A->length, B->coeffs, B->length, n, prec);
    _fmpcb_poly_set_length(Q, n);
    _fmpcb_poly_normalise(Q);
}

