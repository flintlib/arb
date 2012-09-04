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

/* TODO: handle infinities; large shift efficiently */

int
fmprb_contains_fmpr(const fmprb_t x, const fmpr_t y)
{
    if (fmprb_is_exact(x))
    {
        return fmpr_equal(fmprb_midref(x), y);
    }
    else
    {
        fmpr_t t;
        int result = 0;

        fmpr_init(t);

        fmpr_add(t, fmprb_midref(x), fmprb_radref(x),
            FMPR_PREC_EXACT, FMPR_RND_DOWN);

        if (fmpr_cmp(y, t) <= 0)
        {
            fmpr_sub(t, fmprb_midref(x), fmprb_radref(x),
                FMPR_PREC_EXACT, FMPR_RND_DOWN);

            if (fmpr_cmp(y, t) >= 0)
                result = 1;
        }

        fmpr_clear(t);
        return result;
    }
}
