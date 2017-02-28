/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

/* [(a + b) +/- (b - a)] / 2 */

void
arb_set_interval_arf(arb_t x, const arf_t a, const arf_t b, slong prec)
{
    arf_t t;
    int inexact;

    if (arf_is_inf(a) && arf_equal(a, b))
    {
        /* [-inf, -inf] or [+inf, +inf] */
        arf_set(arb_midref(x), a);
        mag_zero(arb_radref(x));
        return;
    }

    arf_init(t);
    arf_sub(t, b, a, MAG_BITS, ARF_RND_UP);

    if (arf_sgn(t) < 0)
    {
        flint_printf("exception: arb_set_interval_arf: endpoints not ordered\n");
        flint_abort();
    }

    arf_get_mag(arb_radref(x), t);

    inexact = arf_add(arb_midref(x), a, b, prec, ARB_RND);
    if (inexact)
        arf_mag_add_ulp(arb_radref(x), arb_radref(x), arb_midref(x), prec);

    arb_mul_2exp_si(x, x, -1);

    arf_clear(t);
}

