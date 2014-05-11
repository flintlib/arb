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

    Copyright (C) 2014 Fredrik Johansson

******************************************************************************/

#include "arb.h"

void
arb_sqrt_ui(arb_t z, ulong x, long prec)
{
    arf_t t;
    arf_init_set_ui(t, x); /* no need to free */
    arb_sqrt_arf(z, t, prec);
}

void
arb_sqrt_fmpz(arb_t z, const fmpz_t x, long prec)
{
    arf_t t;
    arf_init(t);
    arf_set_fmpz(t, x);
    arb_sqrt_arf(z, t, prec);
    arf_clear(t);
}

void
arb_sqrt_arf(arb_t z, const arf_t x, long prec)
{
    if (arf_sgn(x) < 0)
    {
        arb_indeterminate(z);
    }
    else
    {
        int inexact;

        inexact = arf_sqrt(arb_midref(z), x, prec, ARB_RND);

        if (inexact)
            arf_mag_set_ulp(arb_radref(z), arb_midref(z), prec);
        else
            mag_zero(arb_radref(z));
    }
}

void
arb_sqrt(arb_t z, const arb_t x, long prec)
{
    mag_t rx, zr;
    int inexact;

    if (arb_is_exact(x))
    {
        arb_sqrt_arf(z, arb_midref(x), prec);
    }
    else if (arb_contains_negative(x))
    {
        arb_indeterminate(z);
    }
    else
    {
        mag_init(zr);
        mag_init(rx);

        /* rx = upper bound for r / x */
        arf_get_mag_lower(rx, arb_midref(x));
        mag_div(rx, arb_radref(x), rx);

        inexact = arf_sqrt(arb_midref(z), arb_midref(x), prec, ARB_RND);

        /* zr = upper bound for sqrt(x) */
        arf_get_mag(zr, arb_midref(z));
        if (inexact)
            arf_mag_add_ulp(zr, zr, arb_midref(z), prec);

        /* propagated error:   sqrt(x) - sqrt(x-r)
                             = sqrt(x) * [1 - sqrt(1 - r/x)]
                            <= sqrt(x) * 0.5 * (rx + rx^2)  */
        mag_addmul(rx, rx, rx);
        mag_mul(zr, zr, rx);
        mag_mul_2exp_si(zr, zr, -1);

        /* add the rounding error */
        if (inexact)
            arf_mag_add_ulp(arb_radref(z), zr, arb_midref(z), prec);
        else
            mag_swap(arb_radref(z), zr);

        mag_clear(zr);
        mag_clear(rx);
    }
}

