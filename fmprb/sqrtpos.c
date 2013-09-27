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

static __inline__ void
fmprb_nonnegative_part(fmprb_t z, const fmprb_t x, long prec)
{
    if (fmprb_contains_negative(x))
    {
        fmpr_add(fmprb_midref(z), fmprb_midref(x), fmprb_radref(x), prec, FMPR_RND_CEIL);

        if (fmpr_sgn(fmprb_midref(z)) <= 0)
        {
            fmpr_zero(fmprb_radref(z));
        }
        else
        {
            fmpr_mul_2exp_si(fmprb_midref(z), fmprb_midref(z), -1);
            fmpr_set(fmprb_midref(z), fmprb_radref(z));
        }
    }
    else
    {
        fmprb_set(z, x);
    }
}

void
fmprb_sqrtpos(fmprb_t z, const fmprb_t x, long prec)
{
    if (!fmprb_is_finite(x))
    {
        if (fmpr_is_zero(fmprb_radref(x)) && fmpr_is_pos_inf(fmprb_midref(x)))
            fmprb_pos_inf(z);
        else
            fmprb_zero_pm_inf(z);
    }
    else if (fmprb_contains_nonpositive(x))
    {
        fmpr_t t;
        fmpr_init(t);
        fmpr_add(t, fmprb_midref(x), fmprb_radref(x),
            FMPRB_RAD_PREC, FMPR_RND_CEIL);
        if (fmpr_sgn(t) <= 0)
        {
            fmprb_zero(z);
        }
        else
        {
            fmpr_sqrt(t, t, FMPRB_RAD_PREC, FMPR_RND_CEIL);
            fmpr_mul_2exp_si(t, t, -1);
            fmpr_set(fmprb_midref(z), t);
            fmpr_set(fmprb_radref(z), t);
        }
        fmpr_clear(t);
    }
    else
    {
        fmprb_sqrt(z, x, prec);
    }

    fmprb_nonnegative_part(z, z, prec);
}

