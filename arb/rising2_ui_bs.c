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

static void
bsplit(arb_t p, arb_t q, const arb_t x, ulong a, ulong b, long prec)
{
    if (b - a < 8)
    {
        ulong k;
        arb_t t;

        arb_one(p);
        arb_add_ui(q, x, a, prec);

        arb_init(t);

        for (k = a + 1; k < b; k++)
        {
            arb_add_ui(t, x, k, prec);
            arb_mul(p, p, t, prec);
            arb_add(p, p, q, prec);
            arb_mul(q, q, t, prec);
        }

        arb_clear(t);
    }
    else
    {
        arb_t r, s;
        ulong m;

        arb_init(r);
        arb_init(s);

        m = a + (b - a) / 2;
        bsplit(p, q, x, a, m, prec);
        bsplit(r, s, x, m, b, prec);

        arb_mul(p, p, s, prec);
        arb_mul(r, r, q, prec);
        arb_add(p, p, r, prec);
        arb_mul(q, q, s, prec);

        arb_clear(r);
        arb_clear(s);
    }
}

void
arb_rising2_ui_bs(arb_t u, arb_t v, const arb_t x, ulong n, long prec)
{
    if (n == 0)
    {
        arb_zero(v);
        arb_one(u);
    }
    else if (n == 1)
    {
        arb_set(u, x);
        arb_one(v);
    }
    else
    {
        arb_t t;
        long wp = ARF_PREC_ADD(prec, FLINT_BIT_COUNT(n));

        arb_init(t);  /* support aliasing */
        arb_set(t, x);
        bsplit(v, u, t, 0, n, wp);
        arb_clear(t);
    }
}

