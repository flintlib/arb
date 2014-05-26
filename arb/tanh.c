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

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include "arb.h"

void
arb_tanh(arb_t y, const arb_t x, long prec)
{
    arb_t t, u;

    arb_init(t);
    arb_init(u);

    arb_mul_2exp_si(t, x, 1);

    if (arf_sgn(arb_midref(x)) >= 0)
    {
        arb_neg(t, t);
        arb_expm1(t, t, prec + 4);
        arb_add_ui(y, t, 2, prec + 4);
        arb_div(y, t, y, prec);
        arb_neg(y, y);
    }
    else
    {
        arb_expm1(t, t, prec + 4);
        arb_add_ui(y, t, 2, prec + 4);
        arb_div(y, t, y, prec);
    }

    arb_clear(t);
    arb_clear(u);
}

