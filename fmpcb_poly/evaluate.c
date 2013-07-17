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
_fmpcb_poly_evaluate(fmpcb_t res, fmpcb_srcptr f, long len,
                           const fmpcb_t a, long prec)
{
    if (len == 0)
    {
        fmpcb_zero(res);
    }
    else if (len == 1 || fmpcb_is_zero(a))
    {
        fmpcb_set(res, f);
    }
    else
    {
        long i = len - 1;
        fmpcb_t t;

        fmpcb_init(t);
        fmpcb_set(res, f + i);
        for (i = len - 2; i >= 0; i--)
        {
            fmpcb_mul(t, res, a, prec);
            fmpcb_add(res, f + i, t, prec);
        }
        fmpcb_clear(t);
    }
}

void
fmpcb_poly_evaluate(fmpcb_t res, const fmpcb_poly_t f, const fmpcb_t a, long prec)
{
    if (res == a)
    {
        fmpcb_t t;
        fmpcb_init(t);
        _fmpcb_poly_evaluate(t, f->coeffs, f->length, a, prec);
        fmpcb_swap(res, t);
        fmpcb_clear(t);
    }
    else
        _fmpcb_poly_evaluate(res, f->coeffs, f->length, a, prec);
}
