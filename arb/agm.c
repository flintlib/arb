/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_agm(arb_t z, const arb_t x, const arb_t y, slong prec)
{
    arb_t t, u, v, w;

    if (arb_contains_negative(x) || arb_contains_negative(y))
    {
        arb_indeterminate(z);
        return;
    }

    if (arb_is_zero(x) || arb_is_zero(y))
    {
        arb_zero(z);
        return;
    }

    arb_init(t);
    arb_init(u);
    arb_init(v);
    arb_init(w);

    arb_set(t, x);
    arb_set(u, y);

    while (!arb_overlaps(t, u) &&
            !arb_contains_nonpositive(t) &&
            !arb_contains_nonpositive(u))
    {
        arb_add(v, t, u, prec);
        arb_mul_2exp_si(v, v, -1);

        arb_mul(w, t, u, prec);
        arb_sqrt(w, w, prec);

        arb_swap(v, t);
        arb_swap(w, u);
    }

    if (!arb_is_finite(t) || !arb_is_finite(u))
    {
        arb_indeterminate(z);
    }
    else
    {
        arb_union(z, t, u, prec);
    }

    arb_clear(t);
    arb_clear(u);
    arb_clear(v);
    arb_clear(w);
}

