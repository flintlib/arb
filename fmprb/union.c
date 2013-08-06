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

#include "fmprb.h"

void
fmprb_union(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec)
{
    fmpr_t left, right, t;

    if (fmpr_is_nan(fmprb_midref(x)) || fmpr_is_nan(fmprb_midref(y)))
    {
        fmprb_indeterminate(z);
        return;
    }

    if (fmpr_is_pos_inf(fmprb_radref(x)) || fmpr_is_pos_inf(fmprb_radref(y)))
    {
        fmprb_zero_pm_inf(z);
        return;
    }

    fmpr_init(left);
    fmpr_init(right);
    fmpr_init(t);

    fmpr_sub(left, fmprb_midref(x), fmprb_radref(x), prec, FMPR_RND_FLOOR);
    fmpr_sub(t, fmprb_midref(y), fmprb_radref(y), prec, FMPR_RND_FLOOR);
    fmpr_min(left, left, t);

    fmpr_add(right, fmprb_midref(x), fmprb_radref(x), prec, FMPR_RND_CEIL);
    fmpr_add(t, fmprb_midref(y), fmprb_radref(y), prec, FMPR_RND_CEIL);
    fmpr_max(right, right, t);

    fmprb_set_interval_fmpr(z, left, right, prec);

    fmpr_clear(left);
    fmpr_clear(right);
    fmpr_clear(t);
}

