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

#include "arb.h"

void
arb_union(arb_t z, const arb_t x, const arb_t y, slong prec)
{
    arf_t left, right, t, xr, yr;

    if (arf_is_nan(arb_midref(x)) || arf_is_nan(arb_midref(y)))
    {
        arb_indeterminate(z);
        return;
    }

    if (mag_is_inf(arb_radref(x)) || mag_is_inf(arb_radref(y)))
    {
        arb_zero_pm_inf(z);
        return;
    }

    arf_init(left);
    arf_init(right);
    arf_init(t);

    arf_init_set_mag_shallow(xr, arb_radref(x));
    arf_init_set_mag_shallow(yr, arb_radref(y));

    arf_sub(left, arb_midref(x), xr, prec, ARF_RND_FLOOR);
    arf_sub(t, arb_midref(y), yr, prec, ARF_RND_FLOOR);
    arf_min(left, left, t);

    arf_add(right, arb_midref(x), xr, prec, ARF_RND_CEIL);
    arf_add(t, arb_midref(y), yr, prec, ARF_RND_CEIL);
    arf_max(right, right, t);

    arb_set_interval_arf(z, left, right, prec);

    arf_clear(left);
    arf_clear(right);
    arf_clear(t);
}

