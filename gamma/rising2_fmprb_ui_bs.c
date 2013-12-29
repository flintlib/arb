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
bsplit(fmprb_t p, fmprb_t q, const fmprb_t x, ulong a, ulong b, long prec)
{
    if (b - a < 8)
    {
        ulong k;
        fmprb_t t;

        fmprb_one(p);
        fmprb_add_ui(q, x, a, prec);

        fmprb_init(t);

        for (k = a + 1; k < b; k++)
        {
            fmprb_add_ui(t, x, k, prec);
            fmprb_mul(p, p, t, prec);
            fmprb_add(p, p, q, prec);
            fmprb_mul(q, q, t, prec);
        }

        fmprb_clear(t);
    }
    else
    {
        fmprb_t r, s;
        ulong m;

        fmprb_init(r);
        fmprb_init(s);

        m = a + (b - a) / 2;
        bsplit(p, q, x, a, m, prec);
        bsplit(r, s, x, m, b, prec);

        fmprb_mul(p, p, s, prec);
        fmprb_mul(r, r, q, prec);
        fmprb_add(p, p, r, prec);
        fmprb_mul(q, q, s, prec);

        fmprb_clear(r);
        fmprb_clear(s);
    }
}

void
gamma_rising2_fmprb_ui_bs(fmprb_t u, fmprb_t v, const fmprb_t x, ulong n, long prec)
{
    if (n == 0)
    {
        fmprb_zero(v);
        fmprb_one(u);
    }
    else if (n == 1)
    {
        fmprb_set(u, x);
        fmprb_one(v);
    }
    else
    {
        fmprb_t t;
        long wp = FMPR_PREC_ADD(prec, FLINT_BIT_COUNT(n));

        fmprb_init(t);  /* support aliasing */
        fmprb_set(t, x);
        bsplit(v, u, t, 0, n, wp);
        fmprb_clear(t);
    }
}

