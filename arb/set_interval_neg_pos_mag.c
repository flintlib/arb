/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_set_interval_neg_pos_mag(arb_t res, const mag_t a, const mag_t b, slong prec)
{
    arf_t aa, bb;
    int inexact;

    if (mag_is_inf(a) || mag_is_inf(b))
    {
        arb_zero_pm_inf(res);
        return;
    }

    arf_init_set_mag_shallow(aa, a);
    arf_init_set_mag_shallow(bb, b);

    inexact = arf_sub(arb_midref(res), bb, aa, prec, ARB_RND);

    mag_add(arb_radref(res), b, a);

    if (inexact)
        arf_mag_add_ulp(arb_radref(res), arb_radref(res), arb_midref(res), prec);

    arb_mul_2exp_si(res, res, -1);
}

