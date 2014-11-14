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

    Copyright (C) 2014 Fredrik Johansson

******************************************************************************/

#include "arb.h"

void
arb_chebyshev_u_ui(arb_t y, ulong n, const arb_t x, long prec)
{
    arb_t a, b;

    if (n <= 1)
    {
        if (n == 0)
        {
            arb_one(y);
        }
        else
        {
            arb_set_round(y, x, prec);
            arb_mul_2exp_si(y, y, 1);
        }
        return;
    }

    arb_init(a);
    arb_init(b);

    arb_chebyshev_u2_ui(a, b, n / 2, x, prec);

    if (n % 2 == 0)
    {
        arb_add(y, a, b, prec);
        arb_sub(b, a, b, prec);
        arb_mul(y, y, b, prec);
    }
    else
    {
        arb_submul(b, a, x, prec);
        arb_mul(y, a, b, prec);
        arb_mul_2exp_si(y, y, 1);
        arb_neg(y, y);
    }

    arb_clear(a);
    arb_clear(b);
}

