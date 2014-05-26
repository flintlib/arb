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
arb_contains_arf(const arb_t x, const arf_t y)
{
    if (arf_is_nan(y))
    {
        return arf_is_nan(arb_midref(x));
    }
    else if (arf_is_nan(arb_midref(x)))
    {
        return 1;
    }
    else if (arb_is_exact(x))
    {
        return arf_equal(arb_midref(x), y);
    }
    else
    {
        arf_t t;
        arf_struct tmp[3];
        int result;

        arf_init(t);

        /* y >= xm - xr  <=>  0 >= xm - xr - y */
        arf_init_set_shallow(tmp + 0, arb_midref(x));
        arf_init_neg_mag_shallow(tmp + 1,  arb_radref(x));
        arf_init_neg_shallow(tmp + 2, y);

        arf_sum(t, tmp, 3, 30, ARF_RND_DOWN);
        result = (arf_sgn(t) <= 0);

        if (result)
        {
            /* y <= xm + xr  <=>  0 <= xm + xr - y */
            arf_init_set_mag_shallow(tmp + 1,  arb_radref(x));
            arf_sum(t, tmp, 3, 30, ARF_RND_DOWN);
            result = (arf_sgn(t) >= 0);
        }

        arf_clear(t);

        return result;
    }
}

