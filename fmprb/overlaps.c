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
    fmpr_t r, d;
    int result;

    if (fmpr_is_inf(fmprb_radref(x)) || fmpr_is_inf(fmprb_radref(y)))
        return 1;

    fmpr_init(r);
    fmpr_init(d);

    fmpr_sub(d, fmprb_midref(x), fmprb_midref(y), FMPR_PREC_EXACT, FMPR_RND_UP);
    fmpr_abs(d, d);

    fmpr_add(r, fmprb_radref(x), fmprb_radref(y), FMPR_PREC_EXACT, FMPR_RND_UP);

    result = fmpr_cmp(d, r) <= 0;

    fmpr_clear(r);
    fmpr_clear(d);

    return result;
}

