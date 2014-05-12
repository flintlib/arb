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

static __inline__ void
arb_nonnegative_part(arb_t z, const arb_t x, long prec)
{
    if (arb_contains_negative(x))
    {
        arf_t t;
        arf_init(t);

        arf_set_mag(t, arb_radref(x));
        arf_add(arb_midref(z), arb_midref(x), t, MAG_BITS, ARF_RND_CEIL);

        if (arf_sgn(arb_midref(z)) <= 0)
        {
            mag_zero(arb_radref(z));
        }
        else
        {
            arf_mul_2exp_si(arb_midref(z), arb_midref(z), -1);
            arf_get_mag(arb_radref(z), arb_midref(z));

            /* XXX: needed since arf_get_mag is inexact */
            arf_set_mag(arb_midref(z), arb_radref(z));
        }

        arf_clear(t);
    }
    else
    {
        arb_set(z, x);
    }
}

void
arb_sqrtpos(arb_t z, const arb_t x, long prec)
{
    if (!arb_is_finite(x))
    {
        if (mag_is_zero(arb_radref(x)) && arf_is_pos_inf(arb_midref(x)))
            arb_pos_inf(z);
        else
            arb_zero_pm_inf(z);
    }
    else if (arb_contains_nonpositive(x))
    {
        arf_t t;

        arf_init(t);

        arf_set_mag(t, arb_radref(x));
        arf_add(t, arb_midref(x), t, MAG_BITS, ARF_RND_CEIL);

        if (arf_sgn(t) <= 0)
        {
            arb_zero(z);
        }
        else
        {
            arf_sqrt(t, t, MAG_BITS, ARF_RND_CEIL);
            arf_mul_2exp_si(t, t, -1);
            arf_set(arb_midref(z), t);
            arf_get_mag(arb_radref(z), t);
        }

        arf_clear(t);
    }
    else
    {
        arb_sqrt(z, x, prec);
    }

    arb_nonnegative_part(z, z, prec);
}

