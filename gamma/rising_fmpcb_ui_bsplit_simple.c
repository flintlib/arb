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

#include "gamma.h"

/* assumes y and x are not aliased */
static void
bsplit(fmpcb_t y, const fmpcb_t x, ulong a, ulong b, long prec)
{
    if (b - a == 1)
    {
        fmpcb_set_round(y, x, prec);
    }
    else if (b - a <= 10)
    {
        long i;
        fmpcb_t t;
        fmpcb_init(t);

        fmpcb_add_ui(t, x, a, prec);
        fmpcb_add_ui(y, x, a + 1, prec);
        fmpcb_mul(y, y, t, prec);

        for (i = a + 2; i < b; i++)
        {
            fmpcb_add_ui(t, x, i, prec);
            fmpcb_mul(y, y, t, prec);
        }

        fmpcb_clear(t);
    }
    else
    {
        fmpcb_t t, u;
        ulong m = a + (b - a) / 2;

        fmpcb_init(t);
        fmpcb_init(u);

        bsplit(t, x, a, m, prec);
        bsplit(u, x, m, b, prec);

        fmpcb_mul(y, t, u, prec);

        fmpcb_clear(t);
        fmpcb_clear(u);
    }
}

void
gamma_rising_fmpcb_ui_bsplit_simple(fmpcb_t y, const fmpcb_t x, ulong n, long prec)
{
    if (n == 0)
    {
        fmpcb_one(y);
    }
    else if (n == 1)
    {
        fmpcb_set_round(y, x, prec);
    }
    else
    {
        fmpcb_t t;
        long wp = FMPR_PREC_ADD(prec, FLINT_BIT_COUNT(n));

        fmpcb_init(t);
        bsplit(t, x, 0, n, wp);
        fmpcb_set_round(y, t, prec);
        fmpcb_clear(t);
    }
}

