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
arb_chebyshev_t2_ui(arb_t a, arb_t b, ulong n, const arb_t x, long prec)
{
    int i;

    arb_set_round(a, x, prec);
    arb_one(b);

    if (n <= 1)
    {
        if (n == 0)
            arb_swap(a, b);
        return;
    }

    for (i = FLINT_BIT_COUNT(n - 1) - 1; i >= 0; i--)
    {
        if (((n - 1) >> i) & 1)
        {
            arb_mul(b, b, a, prec);
            arb_mul_2exp_si(b, b, 1);
            arb_sub(b, b, x, prec);
            arb_mul(a, a, a, prec);
            arb_mul_2exp_si(a, a, 1);
            arb_sub_ui(a, a, 1, prec);
        }
        else
        {
            arb_mul(a, a, b, prec);
            arb_mul_2exp_si(a, a, 1);
            arb_sub(a, a, x, prec);
            arb_mul(b, b, b, prec);
            arb_mul_2exp_si(b, b, 1);
            arb_sub_ui(b, b, 1, prec);
        }
    }
}

