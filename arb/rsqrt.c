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

    Copyright (C) 2012-2014 Fredrik Johansson

******************************************************************************/

#include <math.h>
#include "arb.h"

void
mag_pow_minus_three_half(mag_t z, const mag_t x)
{
    double t;

    if (mag_is_zero(x))
    {
        mag_inf(z);
    }
    else if (mag_is_inf(x))
    {
        mag_zero(z);
    }
    else
    {
        if (fmpz_is_even(MAG_EXPREF(x)))
        {
            fmpz_mul_si(MAG_EXPREF(z), MAG_EXPREF(x), -3);
            t = MAG_MAN(x) * (1.0 / (LIMB_ONE << MAG_BITS));
        }
        else
        {
            fmpz_add_ui(MAG_EXPREF(z), MAG_EXPREF(x), 1);
            fmpz_mul_si(MAG_EXPREF(z), MAG_EXPREF(z), -3);
            t = MAG_MAN(x) * (0.5 / (LIMB_ONE << MAG_BITS));
        }

        fmpz_tdiv_q_2exp(MAG_EXPREF(z), MAG_EXPREF(z), 1);
        t = 1.0 / (t * sqrt(t));
        t = t * (1.0 + 1e-10);

        mag_set_d_2exp_fmpz(z, t, MAG_EXPREF(z));
    }
}

void
arb_rsqrt(arb_t z, const arb_t x, long prec)
{
    int inexact;

    if (arb_contains_nonpositive(x))
    {
        arb_indeterminate(z);
        return;
    }

    if (arb_is_exact(x))
    {
        inexact = arf_rsqrt(arb_midref(z), arb_midref(x), prec, ARB_RND);

        if (inexact)
            arf_mag_set_ulp(arb_radref(z), arb_midref(z), prec);
        else
            mag_zero(arb_radref(z));
    }
    else
    {
        mag_t t;
        mag_init(t);

        /* error bound: (1/2) (x-r)^(-3/2) * r */
        arb_get_mag_lower(t, x);
        mag_pow_minus_three_half(t, t);
        mag_mul(t, t, arb_radref(x));
        mag_mul_2exp_si(t, t, -1);

        inexact = arf_rsqrt(arb_midref(z), arb_midref(x), prec, ARB_RND);

        if (inexact)
            arf_mag_add_ulp(arb_radref(z), t, arb_midref(z), prec);
        else
            mag_swap(arb_radref(z), t);

        mag_clear(t);
    }
}

void
arb_rsqrt_ui(arb_t z, ulong x, long prec)
{
    arb_t t;
    arb_init(t);
    arb_set_ui(t, x);
    arb_rsqrt(z, t, prec);
    arb_clear(t);
}

