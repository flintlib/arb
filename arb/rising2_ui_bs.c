/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

static void
bsplit(arb_t p, arb_t q, const arb_t x, ulong a, ulong b, slong prec)
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
arb_rising2_ui_bs(arb_t u, arb_t v, const arb_t x, ulong n, slong prec)
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
        slong wp = ARF_PREC_ADD(prec, FLINT_BIT_COUNT(n));

        arb_init(t);  /* support aliasing */
        arb_set(t, x);
        bsplit(v, u, t, 0, n, wp);
        arb_clear(t);
    }
}

