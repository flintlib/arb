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

#include "fmpcb_poly.h"

void
_fmpcb_poly_div_series(fmpcb_struct * Q,
    const fmpcb_struct * A,
    const fmpcb_struct * B, long n, long prec)
{
    fmpcb_struct * Binv = _fmpcb_vec_init(n);

    _fmpcb_poly_inv_series(Binv, B, n, prec);
    _fmpcb_poly_mullow(Q, Binv, n, A, n, n, prec);

    _fmpcb_vec_clear(Binv, n);
}

void
_fmpcb_poly_div(fmpcb_struct * Q,
    const fmpcb_struct * A, long lenA,
    const fmpcb_struct * B, long lenB, long prec)
{
    const long lenQ = lenA - lenB + 1;
    fmpcb_struct * Arev, * Brev;

    Arev = _fmpcb_vec_init(2 * lenQ);
    Brev = Arev + lenQ;

    _fmpcb_poly_reverse(Arev, A + (lenA - lenQ), lenQ, lenQ);

    if (lenB >= lenQ)
    {
        _fmpcb_poly_reverse(Brev, B + (lenB - lenQ), lenQ, lenQ);
    }
    else
    {
        _fmpcb_poly_reverse(Brev, B, lenB, lenB);
        _fmpcb_vec_zero(Brev + lenB, lenQ - lenB);
    }

    _fmpcb_poly_div_series(Q, Arev, Brev, lenQ, prec);
    _fmpcb_poly_reverse(Q, Q, lenQ, lenQ);

    _fmpcb_vec_clear(Arev, 2 * lenQ);
}

void _fmpcb_poly_divrem(fmpcb_struct * Q, fmpcb_struct * R,
    const fmpcb_struct * A, long lenA,
    const fmpcb_struct * B, long lenB, long prec)
{
    const long lenQ = lenA - lenB + 1;
    _fmpcb_poly_div(Q, A, lenA, B, lenB, prec);

    if (lenB > 1)
    {
        if (lenQ >= lenB - 1)
            _fmpcb_poly_mullow(R, Q, lenQ, B, lenB - 1, lenB - 1, prec);
        else
            _fmpcb_poly_mullow(R, B, lenB - 1, Q, lenQ, lenB - 1, prec);
        _fmpcb_vec_sub(R, A, R, lenB - 1, prec);
    }
}

void _fmpcb_poly_rem(fmpcb_struct * R,
    const fmpcb_struct * A, long lenA,
    const fmpcb_struct * B, long lenB, long prec)
{
    const long lenQ = lenA - lenB + 1;
    fmpcb_struct * Q = _fmpcb_vec_init(lenQ);
    _fmpcb_poly_divrem(Q, R, A, lenA, B, lenB, prec);
    _fmpcb_vec_clear(Q, lenQ);
}

void fmpcb_poly_divrem(fmpcb_poly_t Q, fmpcb_poly_t R,
                             const fmpcb_poly_t A, const fmpcb_poly_t B, long prec)
{
    const long lenA = A->length, lenB = B->length;

    if (lenB == 0 || fmpcb_contains_zero(B->coeffs + lenB - 1))
    {
        printf("Exception: division by zero in fmpcb_poly_divrem\n");
        abort();
    }

    if (lenA < lenB)
    {
        fmpcb_poly_set(R, A);
        fmpcb_poly_zero(Q);
        return;
    }

    if (Q == A || Q == B)
    {
        fmpcb_poly_t T;
        fmpcb_poly_init(T);
        fmpcb_poly_divrem(T, R, A, B, prec);
        fmpcb_poly_swap(Q, T);
        fmpcb_poly_clear(T);
        return;
    }

    if (R == A || R == B)
    {
        fmpcb_poly_t U;
        fmpcb_poly_init(U);
        fmpcb_poly_divrem(Q, U, A, B, prec);
        fmpcb_poly_swap(R, U);
        fmpcb_poly_clear(U);
        return;
    }

    fmpcb_poly_fit_length(Q, lenA - lenB + 1);
    fmpcb_poly_fit_length(R, lenB - 1);

    _fmpcb_poly_divrem(Q->coeffs, R->coeffs, A->coeffs, lenA,
                                   B->coeffs, lenB, prec);

    _fmpcb_poly_set_length(Q, lenA - lenB + 1);
    _fmpcb_poly_set_length(R, lenB - 1);
    _fmpcb_poly_normalise(R);
}

