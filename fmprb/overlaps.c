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

#include "fmprb.h"

int
fmprb_overlaps(const fmprb_t x, const fmprb_t y)
{
    fmpr_t t;
    fmpr_struct u[4];
    int result;

    if (fmpr_is_inf(fmprb_radref(x)) || fmpr_is_inf(fmprb_radref(y)))
        return 1;

    if (fmpr_equal(fmprb_midref(x), fmprb_midref(y)))
        return 1;

    fmpr_init(t);
    fmpr_init(u + 0);
    fmpr_init(u + 1);
    fmpr_init(u + 2);
    fmpr_init(u + 3);

    /* |xm - ym| <= xr + yr */

    if (fmpr_cmp(fmprb_midref(x), fmprb_midref(y)) >= 0)
    {
        fmpr_set(u + 0, fmprb_midref(x));
        fmpr_neg(u + 1, fmprb_midref(y));
    }
    else
    {
        fmpr_neg(u + 0, fmprb_midref(x));
        fmpr_set(u + 1, fmprb_midref(y));
    }

    fmpr_neg(u + 2, fmprb_radref(x));
    fmpr_neg(u + 3, fmprb_radref(y));

    fmpr_sum(t, u, 4, 30, FMPR_RND_DOWN);
    result = fmpr_sgn(t) <= 0;

    fmpr_clear(t);
    fmpr_clear(u + 0);
    fmpr_clear(u + 1);
    fmpr_clear(u + 2);
    fmpr_clear(u + 3);

    return result;
}

