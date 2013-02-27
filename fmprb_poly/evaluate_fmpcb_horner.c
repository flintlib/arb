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

#include "fmprb_poly.h"

void
_fmprb_poly_evaluate_fmpcb_horner(fmpcb_t y, const fmprb_struct * f, long len,
                           const fmpcb_t x, long prec)
{
    if (len == 0)
    {
        fmpcb_zero(y);
    }
    else if (len == 1 || fmpcb_is_zero(x))
    {
        fmpcb_set_round_fmprb(y, f, prec);
    }
    else if (len == 2)
    {
        fmpcb_mul_fmprb(y, x, f + 1, prec);
        fmpcb_add_fmprb(y, y, f + 0, prec);
    }
    else
    {
        long i = len - 1;
        fmpcb_t t, u;

        fmpcb_init(t);
        fmpcb_init(u);
        fmpcb_set_fmprb(u, f + i);

        for (i = len - 2; i >= 0; i--)
        {
            fmpcb_mul(t, u, x, prec);
            fmpcb_add_fmprb(u, t, f + i, prec);
        }

        fmpcb_swap(y, u);

        fmpcb_clear(t);
        fmpcb_clear(u);
    }
}

void
fmprb_poly_evaluate_fmpcb_horner(fmpcb_t res, const fmprb_poly_t f, const fmpcb_t a, long prec)
{
    _fmprb_poly_evaluate_fmpcb_horner(res, f->coeffs, f->length, a, prec);
}

