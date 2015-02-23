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

    Copyright (C) 2015 Fredrik Johansson

******************************************************************************/

#include "arb.h"

static void
arb_log1p_tiny(arb_t r, const arb_t z, long prec)
{
    mag_t b, c;
    arb_t t;

    mag_init(b);
    mag_init(c);
    arb_init(t);

    /* if |z| < 1, then |log(1+z) - [z - z^2/2]| <= |z|^3/(1-|z|) */
    arb_get_mag(b, z);
    mag_one(c);
    mag_sub_lower(c, c, b);
    mag_pow_ui(b, b, 3);
    mag_div(b, b, c);

    arb_mul(t, z, z, prec);
    arb_mul_2exp_si(t, t, -1);
    arb_sub(r, z, t, prec);

    if (mag_is_finite(b))
        arb_add_error_mag(r, b);
    else
        arb_indeterminate(r);

    mag_clear(b);
    mag_clear(c);
    arb_clear(t);
}

void
arb_log1p(arb_t r, const arb_t z, long prec)
{
    long magz;

    if (arb_is_zero(z))
    {
        arb_zero(r);
        return;
    }

    magz = arf_abs_bound_lt_2exp_si(arb_midref(z));

    if (magz < -prec)
    {
        arb_log1p_tiny(r, z, prec);
    }
    else
    {
        if (magz < 0)
            arb_add_ui(r, z, 1, prec + (-magz) + 4);
        else
            arb_add_ui(r, z, 1, prec + 4);

        arb_log(r, r, prec);
    }
}

