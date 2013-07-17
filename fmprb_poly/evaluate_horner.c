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
_fmprb_poly_evaluate_horner(fmprb_t y, fmprb_srcptr f, long len,
                           const fmprb_t x, long prec)
{
    if (len == 0)
    {
        fmprb_zero(y);
    }
    else if (len == 1 || fmprb_is_zero(x))
    {
        fmprb_set_round(y, f, prec);
    }
    else if (len == 2)
    {
        fmprb_mul(y, x, f + 1, prec);
        fmprb_add(y, y, f + 0, prec);
    }
    else
    {
        long i = len - 1;
        fmprb_t t, u;

        fmprb_init(t);
        fmprb_init(u);
        fmprb_set(u, f + i);

        for (i = len - 2; i >= 0; i--)
        {
            fmprb_mul(t, u, x, prec);
            fmprb_add(u, f + i, t, prec);
        }

        fmprb_swap(y, u);

        fmprb_clear(t);
        fmprb_clear(u);
    }
}

void
fmprb_poly_evaluate_horner(fmprb_t res, const fmprb_poly_t f, const fmprb_t a, long prec)
{
    _fmprb_poly_evaluate_horner(res, f->coeffs, f->length, a, prec);
}

