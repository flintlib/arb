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

/* assumes y and x are not aliased */
static void
bsplit(arb_t y, const arb_t x, ulong a, ulong b, long prec)
{
    if (b - a == 1)
    {
        arb_set_round(y, x, prec);
    }
    else if (b - a <= 10)
    {
        long i;
        arb_t t;
        arb_init(t);

        arb_add_ui(t, x, a, prec);
        arb_add_ui(y, x, a + 1, prec);
        arb_mul(y, y, t, prec);

        for (i = a + 2; i < b; i++)
        {
            arb_add_ui(t, x, i, prec);
            arb_mul(y, y, t, prec);
        }

        arb_clear(t);
    }
    else
    {
        arb_t t, u;
        ulong m = a + (b - a) / 2;

        arb_init(t);
        arb_init(u);

        bsplit(t, x, a, m, prec);
        bsplit(u, x, m, b, prec);

        arb_mul(y, t, u, prec);

        arb_clear(t);
        arb_clear(u);
    }
}

void
arb_rising_ui_bs(arb_t y, const arb_t x, ulong n, long prec)
{
    if (n == 0)
    {
        arb_one(y);
    }
    else if (n == 1)
    {
        arb_set_round(y, x, prec);
    }
    else
    {
        arb_t t;
        long wp = ARF_PREC_ADD(prec, FLINT_BIT_COUNT(n));

        arb_init(t);
        bsplit(t, x, 0, n, wp);
        arb_set_round(y, t, prec);
        arb_clear(t);
    }
}

