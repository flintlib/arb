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
fmprb_contains_fmpr(const fmprb_t x, const fmpr_t y)
{
    if (fmpr_is_nan(y))
    {
        return fmpr_is_nan(fmprb_midref(x));
    }
    else if (fmpr_is_nan(fmprb_midref(x)))
    {
        return 1;
    }
    else if (fmprb_is_exact(x))
    {
        return fmpr_equal(fmprb_midref(x), y);
    }
    else
    {
        fmpr_t t;
        fmpr_struct tmp[3];
        int result;

        fmpr_init(t);
        fmpr_init(tmp + 0);
        fmpr_init(tmp + 1);
        fmpr_init(tmp + 2);

        /* y >= xm - xr  <=>  0 >= xm - xr - y */
        fmpr_set(tmp + 0, fmprb_midref(x));
        fmpr_neg(tmp + 1, fmprb_radref(x));
        fmpr_neg(tmp + 2, y);
        fmpr_sum(t, tmp, 3, 30, FMPR_RND_DOWN);
        result = (fmpr_sgn(t) <= 0);

        if (result)
        {
            /* y <= xm + xr  <=>  0 <= xm + xr - y */
            fmpr_neg(tmp + 1, tmp + 1);
            fmpr_sum(t, tmp, 3, 30, FMPR_RND_DOWN);
            result = (fmpr_sgn(t) >= 0);
        }

        fmpr_clear(t);
        fmpr_clear(tmp + 0);
        fmpr_clear(tmp + 1);
        fmpr_clear(tmp + 2);

        return result;
    }
}

