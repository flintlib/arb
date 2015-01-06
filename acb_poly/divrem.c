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

    Copyright (C) 2011 William Hart
    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "acb_poly.h"

/* TODO: tighten this code */

void
_acb_poly_div(acb_ptr Q,
    acb_srcptr A, long lenA,
    acb_srcptr B, long lenB, long prec)
{
    long lenQ, lenB2;
    acb_ptr Arev, Brev;

    lenQ = lenA - lenB + 1;

    Arev = _acb_vec_init(2 * lenQ);
    Brev = Arev + lenQ;

    _acb_poly_reverse(Arev, A + (lenA - lenQ), lenQ, lenQ);

    if (lenB >= lenQ)
    {
        _acb_poly_reverse(Brev, B + (lenB - lenQ), lenQ, lenQ);
        lenB2 = lenQ;
    }
    else
    {
        _acb_poly_reverse(Brev, B, lenB, lenB);
        lenB2 = lenB;
    }

    _acb_poly_div_series(Q, Arev, lenQ, Brev, lenB2, lenQ, prec);
    _acb_poly_reverse(Q, Q, lenQ, lenQ);

    _acb_vec_clear(Arev, 2 * lenQ);
}

void _acb_poly_divrem(acb_ptr Q, acb_ptr R,
    acb_srcptr A, long lenA,
    acb_srcptr B, long lenB, long prec)
{
    const long lenQ = lenA - lenB + 1;
    _acb_poly_div(Q, A, lenA, B, lenB, prec);

    if (lenB > 1)
    {
        if (lenQ >= lenB - 1)
            _acb_poly_mullow(R, Q, lenQ, B, lenB - 1, lenB - 1, prec);
        else
            _acb_poly_mullow(R, B, lenB - 1, Q, lenQ, lenB - 1, prec);
        _acb_vec_sub(R, A, R, lenB - 1, prec);
    }
}

void _acb_poly_rem(acb_ptr R,
    acb_srcptr A, long lenA,
    acb_srcptr B, long lenB, long prec)
{
    const long lenQ = lenA - lenB + 1;
    acb_ptr Q = _acb_vec_init(lenQ);
    _acb_poly_divrem(Q, R, A, lenA, B, lenB, prec);
    _acb_vec_clear(Q, lenQ);
}

int acb_poly_divrem(acb_poly_t Q, acb_poly_t R,
                             const acb_poly_t A, const acb_poly_t B, long prec)
{
    const long lenA = A->length, lenB = B->length;

    if (lenB == 0 || acb_contains_zero(B->coeffs + lenB - 1))
    {
        return 0;
    }

    if (lenA < lenB)
    {
        acb_poly_set(R, A);
        acb_poly_zero(Q);
        return 1;
    }

    if (Q == A || Q == B)
    {
        acb_poly_t T;
        acb_poly_init(T);
        acb_poly_divrem(T, R, A, B, prec);
        acb_poly_swap(Q, T);
        acb_poly_clear(T);
        return 1;
    }

    if (R == A || R == B)
    {
        acb_poly_t U;
        acb_poly_init(U);
        acb_poly_divrem(Q, U, A, B, prec);
        acb_poly_swap(R, U);
        acb_poly_clear(U);
        return 1;
    }

    acb_poly_fit_length(Q, lenA - lenB + 1);
    acb_poly_fit_length(R, lenB - 1);

    _acb_poly_divrem(Q->coeffs, R->coeffs, A->coeffs, lenA,
                                   B->coeffs, lenB, prec);

    _acb_poly_set_length(Q, lenA - lenB + 1);
    _acb_poly_set_length(R, lenB - 1);
    _acb_poly_normalise(R);

    return 1;
}
