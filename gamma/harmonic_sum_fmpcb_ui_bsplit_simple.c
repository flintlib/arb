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

static void
bsplit(fmpcb_t p, fmpcb_t q, const fmpcb_t x, ulong a, ulong b, long prec)
{
    if (b - a < 8)
    {
        ulong k;
        fmpcb_t t;

        fmpcb_one(p);
        fmpcb_add_ui(q, x, a, prec);

        fmpcb_init(t);

        for (k = a + 1; k < b; k++)
        {
            fmpcb_add_ui(t, x, k, prec);
            fmpcb_mul(p, p, t, prec);
            fmpcb_add(p, p, q, prec);
            fmpcb_mul(q, q, t, prec);
        }

        fmpcb_clear(t);
    }
    else
    {
        fmpcb_t r, s;
        ulong m;

        fmpcb_init(r);
        fmpcb_init(s);

        m = a + (b - a) / 2;
        bsplit(p, q, x, a, m, prec);
        bsplit(r, s, x, m, b, prec);

        fmpcb_mul(p, p, s, prec);
        fmpcb_mul(r, r, q, prec);
        fmpcb_add(p, p, r, prec);
        fmpcb_mul(q, q, s, prec);

        fmpcb_clear(r);
        fmpcb_clear(s);
    }
}

void
gamma_harmonic_sum_fmpcb_ui_bsplit_simple(fmpcb_t y, const fmpcb_t x, ulong n, long prec)
{
    if (n == 0)
    {
        fmpcb_zero(y);
    }
    else if (n == 1)
    {
        fmpcb_inv(y, x, prec);
    }
    else
    {
        fmpcb_t p, q;
        long wp;

        wp = FMPR_PREC_ADD(prec, FLINT_BIT_COUNT(n));

        fmpcb_init(p);
        fmpcb_init(q);

        bsplit(p, q, x, 0, n, wp);
        fmpcb_div(y, p, q, prec);

        fmpcb_clear(p);
        fmpcb_clear(q);
    }
}

