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
fmprb_contains(const fmprb_t x, const fmprb_t y)
{
    fmpr_t t;
    fmpr_t u;
    int left_ok, right_ok;

    fmpr_init(t);
    fmpr_init(u);

    fmpr_sub(t, fmprb_midref(x), fmprb_radref(x), 30, FMPR_RND_CEIL);
    fmpr_sub(u, fmprb_midref(y), fmprb_radref(y), 30, FMPR_RND_FLOOR);
    left_ok = fmpr_cmp(t, u) <= 0;

    if (!left_ok)
    {
        fmpr_sub(t, fmprb_midref(x), fmprb_radref(x), FMPR_PREC_EXACT, FMPR_RND_CEIL);
        fmpr_sub(u, fmprb_midref(y), fmprb_radref(y), FMPR_PREC_EXACT, FMPR_RND_FLOOR);
        left_ok = fmpr_cmp(t, u) <= 0;
    }

    fmpr_add(t, fmprb_midref(x), fmprb_radref(x), 30, FMPR_RND_FLOOR);
    fmpr_add(u, fmprb_midref(y), fmprb_radref(y), 30, FMPR_RND_CEIL);
    right_ok = (fmpr_cmp(t, u) >= 0);

    if (!right_ok)
    {
        fmpr_add(t, fmprb_midref(x), fmprb_radref(x), FMPR_PREC_EXACT, FMPR_RND_FLOOR);
        fmpr_add(u, fmprb_midref(y), fmprb_radref(y), FMPR_PREC_EXACT, FMPR_RND_CEIL);
        right_ok = (fmpr_cmp(t, u) >= 0);
    }

    fmpr_clear(t);
    fmpr_clear(u);

    return left_ok && right_ok;
}

