/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb.h"

static void
bsplit(acb_t p, acb_t q, const acb_t x, ulong a, ulong b, slong prec)
{
    if (b - a < 8)
    {
        ulong k;
        acb_t t;

        acb_one(p);
        acb_add_ui(q, x, a, prec);

        acb_init(t);

        for (k = a + 1; k < b; k++)
        {
            acb_add_ui(t, x, k, prec);
            acb_mul(p, p, t, prec);
            acb_add(p, p, q, prec);
            acb_mul(q, q, t, prec);
        }

        acb_clear(t);
    }
    else
    {
        acb_t r, s;
        ulong m;

        acb_init(r);
        acb_init(s);

        m = a + (b - a) / 2;
        bsplit(p, q, x, a, m, prec);
        bsplit(r, s, x, m, b, prec);

        acb_mul(p, p, s, prec);
        acb_mul(r, r, q, prec);
        acb_add(p, p, r, prec);
        acb_mul(q, q, s, prec);

        acb_clear(r);
        acb_clear(s);
    }
}

void
acb_rising2_ui_bs(acb_t u, acb_t v, const acb_t x, ulong n, slong prec)
{
    if (n == 0)
    {
        acb_zero(v);
        acb_one(u);
    }
    else if (n == 1)
    {
        acb_set(u, x);
        acb_one(v);
    }
    else
    {
        acb_t t;
        slong wp = ARF_PREC_ADD(prec, FLINT_BIT_COUNT(n));

        acb_init(t);  /* support aliasing */
        acb_set(t, x);
        bsplit(v, u, t, 0, n, wp);
        acb_clear(t);
    }
}

