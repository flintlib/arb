/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb.h"

/* assumes y and x are not aliased */
static void
bsplit(acb_t y, const acb_t x, ulong a, ulong b, slong prec)
{
    if (b - a == 1)
    {
        acb_set_round(y, x, prec);
    }
    else if (b - a <= 10)
    {
        slong i;
        acb_t t;
        acb_init(t);

        acb_add_ui(t, x, a, prec);
        acb_add_ui(y, x, a + 1, prec);
        acb_mul(y, y, t, prec);

        for (i = a + 2; i < b; i++)
        {
            acb_add_ui(t, x, i, prec);
            acb_mul(y, y, t, prec);
        }

        acb_clear(t);
    }
    else
    {
        acb_t t, u;
        ulong m = a + (b - a) / 2;

        acb_init(t);
        acb_init(u);

        bsplit(t, x, a, m, prec);
        bsplit(u, x, m, b, prec);

        acb_mul(y, t, u, prec);

        acb_clear(t);
        acb_clear(u);
    }
}

void
acb_rising_ui_bs(acb_t y, const acb_t x, ulong n, slong prec)
{
    if (n == 0)
    {
        acb_one(y);
    }
    else if (n == 1)
    {
        acb_set_round(y, x, prec);
    }
    else
    {
        acb_t t;
        slong wp = ARF_PREC_ADD(prec, FLINT_BIT_COUNT(n));

        acb_init(t);
        bsplit(t, x, 0, n, wp);
        acb_set_round(y, t, prec);
        acb_clear(t);
    }
}

