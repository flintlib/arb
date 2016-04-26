/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

int
arb_overlaps(const arb_t x, const arb_t y)
{
    arf_t t;
    arf_struct u[4];
    int result;

    if (mag_is_inf(arb_radref(x)) || mag_is_inf(arb_radref(y)))
        return 1;

    if (arf_equal(arb_midref(x), arb_midref(y)))
        return 1;

    arf_init(t);

    /* |xm - ym| <= xr + yr */

    if (arf_cmp(arb_midref(x), arb_midref(y)) >= 0)
    {
        arf_init_set_shallow(u + 0, arb_midref(x));
        arf_init_neg_shallow(u + 1, arb_midref(y));
    }
    else
    {
        arf_init_neg_shallow(u + 0, arb_midref(x));
        arf_init_set_shallow(u + 1, arb_midref(y));
    }

    arf_init_neg_mag_shallow(u + 2, arb_radref(x));
    arf_init_neg_mag_shallow(u + 3, arb_radref(y));

    arf_sum(t, u, 4, MAG_BITS, ARF_RND_DOWN);
    result = arf_sgn(t) <= 0;

    arf_clear(t);

    return result;
}

