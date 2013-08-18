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

#define TRIM_PADDING 16

void
fmprb_trim(fmprb_t y, const fmprb_t x)
{
    if (fmpr_is_zero(fmprb_radref(x)) || fmpr_is_special(fmprb_midref(x)))
    {
        fmprb_set(y, x);
    }
    else if (fmpr_is_special(fmprb_radref(x)))
    {
        /* midpoint must be finite, so set to 0 +/- inf */
        fmprb_zero_pm_inf(y);
    }
    else
    {
        long bits, accuracy;

        bits = fmprb_bits(x);
        accuracy = fmprb_rel_accuracy_bits(x);

        if (accuracy < -TRIM_PADDING)
        {
            /* set to 0 +/- rad */
            if (fmpr_sgn(fmprb_midref(x)) < 0)
                fmpr_sub(fmprb_radref(y), fmprb_radref(x), fmprb_midref(x), FMPRB_RAD_PREC, FMPR_RND_UP);
            else
                fmpr_add(fmprb_radref(y), fmprb_radref(x), fmprb_midref(x), FMPRB_RAD_PREC, FMPR_RND_UP);
            fmpr_zero(fmprb_midref(y));
        }
        else if (accuracy < bits - TRIM_PADDING)
        {
            fmprb_set_round(y, x, FLINT_MAX(0, accuracy) + TRIM_PADDING);
        }
        else
        {
            fmprb_set(y, x);
        }
    }
}

