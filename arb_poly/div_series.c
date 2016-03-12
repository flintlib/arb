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

void 
_arb_poly_div_series(arb_ptr Q, arb_srcptr A, slong Alen,
    arb_srcptr B, slong Blen, slong n, slong prec)
{
    Alen = FLINT_MIN(Alen, n);
    Blen = FLINT_MIN(Blen, n);

    if (Blen == 1)
    {
        _arb_vec_scalar_div(Q, A, Alen, B, prec);
        _arb_vec_zero(Q + Alen, n - Alen);
    }
    else if (n == 2)
    {
        if (Alen == 1)
        {
            arb_div(Q, A, B, prec);
            arb_div(Q + 1, Q, B, prec);
            arb_mul(Q + 1, Q + 1, B + 1, prec);
            arb_neg(Q + 1, Q + 1);
        }
        else
        {
            arb_div(Q, A, B, prec);
            arb_mul(Q + 1, Q, B + 1, prec);
            arb_sub(Q + 1, A + 1, Q + 1, prec);
            arb_div(Q + 1, Q + 1, B, prec);
        }
    }
    else if (Blen == 2 || n <= 10)
    {
        /* The basecase algorithm is faster for much larger Blen and n than
           this, but unfortunately has worse numerical stability. */
        slong i, j;
        arb_t q;

        arb_init(q);

        arb_inv(q, B, prec);
        arb_div(Q, A, B, prec);

        for (i = 1; i < n; i++)
        {
            arb_mul(Q + i, B + 1, Q + i - 1, prec);

            for (j = 2; j < FLINT_MIN(i + 1, Blen); j++)
                arb_addmul(Q + i, B + j, Q + i - j, prec);

            if (i < Alen)
                arb_sub(Q + i, A + i, Q + i, prec);
            else
                arb_neg(Q + i, Q + i);

            if (!arb_is_one(q))
                arb_mul(Q + i, Q + i, q, prec);
        }

        arb_clear(q);
    }
    else
    {
        arb_ptr Binv;
        Binv = _arb_vec_init(n);
        _arb_poly_inv_series(Binv, B, Blen, n, prec);
        _arb_poly_mullow(Q, Binv, n, A, Alen, n, prec);
        _arb_vec_clear(Binv, n);
    }
}

void
arb_poly_div_series(arb_poly_t Q, const arb_poly_t A, const arb_poly_t B, slong n, slong prec)
{
    if (n == 0)
    {
        arb_poly_zero(Q);
        return;
    }

    if (B->length == 0)
    {
        arb_poly_fit_length(Q, n);
        _arb_vec_indeterminate(Q->coeffs, n);
        _arb_poly_set_length(Q, n);
        return;
    }

    if (A->length == 0)
    {
        arb_poly_zero(Q);
        return;
    }

    if (Q == A || Q == B)
    {
        arb_poly_t t;
        arb_poly_init(t);
        arb_poly_div_series(t, A, B, n, prec);
        arb_poly_swap(Q, t);
        arb_poly_clear(t);
        return;
    }

    arb_poly_fit_length(Q, n);
    _arb_poly_div_series(Q->coeffs, A->coeffs, A->length, B->coeffs, B->length, n, prec);
    _arb_poly_set_length(Q, n);
    _arb_poly_normalise(Q);
}

