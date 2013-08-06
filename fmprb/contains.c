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

    Copyright (C) 2012, 2013 Fredrik Johansson

******************************************************************************/

#include "fmprb.h"

int
fmprb_contains(const fmprb_t x, const fmprb_t y)
{
    fmpr_t t;
    fmpr_t u;
    fmpr_struct tmp[4];
    int left_ok, right_ok;

    if (fmprb_is_exact(y))
        return fmprb_contains_fmpr(x, fmprb_midref(y));

    if (fmpr_is_nan(fmprb_midref(y)))
        return fmpr_is_nan(fmprb_midref(x));

    fmpr_init(t);
    fmpr_init(u);

    /* fast check */
    fmpr_sub(t, fmprb_midref(x), fmprb_radref(x), 30, FMPR_RND_CEIL);
    fmpr_sub(u, fmprb_midref(y), fmprb_radref(y), 30, FMPR_RND_FLOOR);
    left_ok = fmpr_cmp(t, u) <= 0;

    /* exact check */
    if (!left_ok)
    {
        fmpr_init(tmp + 0);
        fmpr_init(tmp + 1);
        fmpr_init(tmp + 2);
        fmpr_init(tmp + 3);

        fmpr_set(tmp + 0, fmprb_midref(x));
        fmpr_neg(tmp + 1, fmprb_radref(x));
        fmpr_neg(tmp + 2, fmprb_midref(y));
        fmpr_set(tmp + 3, fmprb_radref(y));

        fmpr_sum(t, tmp, 4, 30, FMPR_RND_DOWN);
        left_ok = fmpr_sgn(t) <= 0;

        fmpr_clear(tmp + 0);
        fmpr_clear(tmp + 1);
        fmpr_clear(tmp + 2);
        fmpr_clear(tmp + 3);
    }

    /* fast check */
    fmpr_add(t, fmprb_midref(x), fmprb_radref(x), 30, FMPR_RND_FLOOR);
    fmpr_add(u, fmprb_midref(y), fmprb_radref(y), 30, FMPR_RND_CEIL);
    right_ok = (fmpr_cmp(t, u) >= 0);

    /* exact check */
    if (!right_ok)
    {
        fmpr_init(tmp + 0);
        fmpr_init(tmp + 1);
        fmpr_init(tmp + 2);
        fmpr_init(tmp + 3);

        fmpr_set(tmp + 0, fmprb_midref(x));
        fmpr_set(tmp + 1, fmprb_radref(x));
        fmpr_neg(tmp + 2, fmprb_midref(y));
        fmpr_neg(tmp + 3, fmprb_radref(y));

        fmpr_sum(t, tmp, 4, 30, FMPR_RND_DOWN);
        right_ok = fmpr_sgn(t) >= 0;

        fmpr_clear(tmp + 0);
        fmpr_clear(tmp + 1);
        fmpr_clear(tmp + 2);
        fmpr_clear(tmp + 3);
    }

    fmpr_clear(t);
    fmpr_clear(u);

    return left_ok && right_ok;
}

