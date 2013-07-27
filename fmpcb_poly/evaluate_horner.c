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

    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "fmpcb_poly.h"

void
_fmpcb_poly_evaluate_horner(fmpcb_t y, fmpcb_srcptr f, long len,
                           const fmpcb_t x, long prec)
{
    if (len == 0)
    {
        fmpcb_zero(y);
    }
    else if (len == 1 || fmpcb_is_zero(x))
    {
        fmpcb_set_round(y, f, prec);
    }
    else if (len == 2)
    {
        fmpcb_mul(y, x, f + 1, prec);
        fmpcb_add(y, y, f + 0, prec);
    }
    else
    {
        long i = len - 1;
        fmpcb_t t, u;

        fmpcb_init(t);
        fmpcb_init(u);
        fmpcb_set(u, f + i);

        for (i = len - 2; i >= 0; i--)
        {
            fmpcb_mul(t, u, x, prec);
            fmpcb_add(u, f + i, t, prec);
        }

        fmpcb_swap(y, u);

        fmpcb_clear(t);
        fmpcb_clear(u);
    }
}

void
fmpcb_poly_evaluate_horner(fmpcb_t res, const fmpcb_poly_t f, const fmpcb_t a, long prec)
{
    _fmpcb_poly_evaluate_horner(res, f->coeffs, f->length, a, prec);
}

