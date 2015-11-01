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

#include "acb.h"

void
acb_chebyshev_t2_ui(acb_t a, acb_t b, ulong n, const acb_t x, long prec)
{
    int i;

    acb_set_round(a, x, prec);
    acb_one(b);

    if (n <= 1)
    {
        if (n == 0)
            acb_swap(a, b);
        return;
    }

    for (i = FLINT_BIT_COUNT(n - 1) - 1; i >= 0; i--)
    {
        if (((n - 1) >> i) & 1)
        {
            acb_mul(b, b, a, prec);
            acb_mul_2exp_si(b, b, 1);
            acb_sub(b, b, x, prec);
            acb_mul(a, a, a, prec);
            acb_mul_2exp_si(a, a, 1);
            acb_sub_ui(a, a, 1, prec);
        }
        else
        {
            acb_mul(a, a, b, prec);
            acb_mul_2exp_si(a, a, 1);
            acb_sub(a, a, x, prec);
            acb_mul(b, b, b, prec);
            acb_mul_2exp_si(b, b, 1);
            acb_sub_ui(b, b, 1, prec);
        }
    }
}

