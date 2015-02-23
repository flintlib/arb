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

#include "acb.h"

static void
acb_log1p_tiny(acb_t r, const acb_t z, long prec)
{
    mag_t b, c;
    acb_t t;
    int real;

    mag_init(b);
    mag_init(c);
    acb_init(t);

    real = acb_is_real(z);

    /* if |z| < 1, then |log(1+z) - [z - z^2/2]| <= |z|^3/(1-|z|) */
    acb_get_mag(b, z);
    mag_one(c);
    mag_sub_lower(c, c, b);
    mag_pow_ui(b, b, 3);
    mag_div(b, b, c);

    acb_mul(t, z, z, prec);
    acb_mul_2exp_si(t, t, -1);
    acb_sub(r, z, t, prec);

    if (real && mag_is_finite(b))
        arb_add_error_mag(acb_realref(r), b);
    else
        acb_add_error_mag(r, b);

    mag_clear(b);
    mag_clear(c);
    acb_clear(t);
}

void
acb_log1p(acb_t r, const acb_t z, long prec)
{
    long magz, magx, magy;

    if (acb_is_zero(z))
    {
        acb_zero(r);
        return;
    }

    magx = arf_abs_bound_lt_2exp_si(arb_midref(acb_realref(z)));
    magy = arf_abs_bound_lt_2exp_si(arb_midref(acb_imagref(z)));
    magz = FLINT_MAX(magx, magy);

    if (magz < -prec)
    {
        acb_log1p_tiny(r, z, prec);
    }
    else
    {
        if (magz < 0)
            acb_add_ui(r, z, 1, prec + (-magz) + 4);
        else
            acb_add_ui(r, z, 1, prec + 4);

        acb_log(r, r, prec);
    }
}

