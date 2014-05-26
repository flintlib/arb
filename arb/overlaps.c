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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

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

    arf_sum(t, u, 4, 30, ARF_RND_DOWN);
    result = arf_sgn(t) <= 0;

    arf_clear(t);

    return result;
}

