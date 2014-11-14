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
arb_chebyshev_t_ui(arb_t y, ulong n, const arb_t x, long prec)
{
    int i, r;

    if (n <= 1)
    {
        if (n == 0)
            arb_one(y);
        else
            arb_set_round(y, x, prec);
        return;
    }

    count_trailing_zeros(r, n);

    if ((n >> r) == 1)
    {
        arb_mul(y, x, x, prec);
        arb_mul_2exp_si(y, y, 1);
        arb_sub_ui(y, y, 1, prec);
        r -= 1;
    }
    else
    {
        /* we only need one value, so break out final iteration */
        arb_t t, u;
        arb_init(t);
        arb_init(u);
        arb_chebyshev_t2_ui(t, u, (n >> (r + 1)) + 1, x, prec);
        arb_mul(t, t, u, prec);
        arb_mul_2exp_si(t, t, 1);
        arb_sub(y, t, x, prec);
        arb_clear(t);
        arb_clear(u);
    }

    for (i = 0; i < r; i++)
    {
        arb_mul(y, y, y, prec);
        arb_mul_2exp_si(y, y, 1);
        arb_sub_ui(y, y, 1, prec);
    }
}

