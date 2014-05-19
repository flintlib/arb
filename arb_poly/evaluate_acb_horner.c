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

#include "arb_poly.h"

void
_arb_poly_evaluate_acb_horner(acb_t y, arb_srcptr f, long len,
                           const acb_t x, long prec)
{
    if (len == 0)
    {
        acb_zero(y);
    }
    else if (len == 1 || acb_is_zero(x))
    {
        acb_set_round_arb(y, f, prec);
    }
    else if (len == 2)
    {
        acb_mul_arb(y, x, f + 1, prec);
        acb_add_arb(y, y, f + 0, prec);
    }
    else
    {
        long i = len - 1;
        acb_t t, u;

        acb_init(t);
        acb_init(u);
        acb_set_arb(u, f + i);

        for (i = len - 2; i >= 0; i--)
        {
            acb_mul(t, u, x, prec);
            acb_add_arb(u, t, f + i, prec);
        }

        acb_swap(y, u);

        acb_clear(t);
        acb_clear(u);
    }
}

void
arb_poly_evaluate_acb_horner(acb_t res, const arb_poly_t f, const acb_t a, long prec)
{
    _arb_poly_evaluate_acb_horner(res, f->coeffs, f->length, a, prec);
}

