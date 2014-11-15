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
arb_div_arf(arb_t z, const arb_t x, const arf_t y, long prec)
{
    mag_t zr, ym;
    int inexact;

    if (arf_is_zero(y))
    {
        if (arf_is_nan(arb_midref(x)))
            arb_indeterminate(z);
        else
            arb_zero_pm_inf(z);
    }
    else if (arb_is_exact(x))
    {
        inexact = arf_div(arb_midref(z), arb_midref(x), y, prec, ARB_RND);

        if (inexact)
            arf_mag_set_ulp(arb_radref(z), arb_midref(z), prec);
        else
            mag_zero(arb_radref(z));
    }
    else if (mag_is_inf(arb_radref(x)))
    {
        arf_div(arb_midref(z), arb_midref(x), y, prec, ARB_RND);
        mag_inf(arb_radref(z));
    }
    else
    {
        mag_init(ym);
        mag_init(zr);

        arf_get_mag_lower(ym, y);
        mag_div(zr, arb_radref(x), ym);

        inexact = arf_div(arb_midref(z), arb_midref(x), y, prec, ARB_RND);

        if (inexact)
            arf_mag_add_ulp(arb_radref(z), zr, arb_midref(z), prec);
        else
            mag_swap(arb_radref(z), zr);

        mag_clear(ym);
        mag_clear(zr);
    }
}

void
arb_div(arb_t z, const arb_t x, const arb_t y, long prec)
{
    mag_t zr, xm, ym, yl, yw;
    int inexact;

    if (arb_is_exact(y))
    {
        arb_div_arf(z, x, arb_midref(y), prec);
    }
    else if (mag_is_inf(arb_radref(x)) || mag_is_inf(arb_radref(y)))
    {
        arf_div(arb_midref(z), arb_midref(x), arb_midref(y), prec, ARB_RND);
        mag_inf(arb_radref(z));
    }
    else
    {
        mag_init_set_arf(xm, arb_midref(x));
        mag_init_set_arf(ym, arb_midref(y));
        mag_init(zr);
        mag_init(yl);
        mag_init(yw);

        /* (|x|*yrad + |y|*xrad)/(y*(|y|-yrad)) */
        mag_mul(zr, xm, arb_radref(y));
        mag_addmul(zr, ym, arb_radref(x));
        arb_get_mag_lower(yw, y);

        arf_get_mag_lower(yl, arb_midref(y));
        mag_mul_lower(yl, yl, yw);

        mag_div(zr, zr, yl);

        inexact = arf_div(arb_midref(z), arb_midref(x), arb_midref(y), prec, ARB_RND);

        if (inexact)
            arf_mag_add_ulp(arb_radref(z), zr, arb_midref(z), prec);
        else
            mag_swap(arb_radref(z), zr);

        mag_clear(xm);
        mag_clear(ym);
        mag_clear(zr);
        mag_clear(yl);
        mag_clear(yw);
    }
}

void
arb_div_si(arb_t z, const arb_t x, long y, long prec)
{
    arf_t t;
    arf_init_set_si(t, y); /* no need to free */
    arb_div_arf(z, x, t, prec);
}

void
arb_div_ui(arb_t z, const arb_t x, ulong y, long prec)
{
    arf_t t;
    arf_init_set_ui(t, y); /* no need to free */
    arb_div_arf(z, x, t, prec);
}

void
arb_div_fmpz(arb_t z, const arb_t x, const fmpz_t y, long prec)
{
    arf_t t;

    if (!COEFF_IS_MPZ(*y))
    {
        arf_init_set_si(t, *y); /* no need to free */
        arb_div_arf(z, x, t, prec);
    }
    else
    {
        arf_init(t);
        arf_set_fmpz(t, y);
        arb_div_arf(z, x, t, prec);
        arf_clear(t);
    }
}


void
arb_fmpz_div_fmpz(arb_t z, const fmpz_t x, const fmpz_t y, long prec)
{
    int inexact;

    inexact = arf_fmpz_div_fmpz(arb_midref(z), x, y, prec, ARB_RND);

    if (inexact)
        arf_mag_set_ulp(arb_radref(z), arb_midref(z), prec);
    else
        mag_zero(arb_radref(z));
}

void
arb_ui_div(arb_t z, ulong x, const arb_t y, long prec)
{
    arb_t t;
    arb_init(t);
    arb_set_ui(t, x);
    arb_div(z, t, y, prec);
    arb_clear(t);
}

