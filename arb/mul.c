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
arb_mul_arf(arb_t z, const arb_t x, const arf_t y, long prec)
{
    mag_t zr, ym;
    int inexact;

    if (arb_is_exact(x))
    {
        inexact = arf_mul(arb_midref(z), arb_midref(x), y, prec, ARB_RND);

        if (inexact)
            arf_mag_set_ulp(arb_radref(z), arb_midref(z), prec);
        else
            mag_zero(arb_radref(z));
    }
    else if (ARB_IS_LAGOM(x) && ARF_IS_LAGOM(y) && ARB_IS_LAGOM(z))
    {
        mag_fast_init_set_arf(ym, y);
        mag_init(zr);

        mag_fast_mul(zr, ym, arb_radref(x));

        inexact = arf_mul(arb_midref(z), arb_midref(x), y, prec, ARB_RND);

        if (inexact)
            arf_mag_fast_add_ulp(zr, zr, arb_midref(z), prec);

        *arb_radref(z) = *zr;
    }
    else
    {
        mag_init_set_arf(ym, y);
        mag_init(zr);

        mag_mul(zr, ym, arb_radref(x));

        inexact = arf_mul(arb_midref(z), arb_midref(x), y, prec, ARB_RND);

        if (inexact)
            arf_mag_add_ulp(arb_radref(z), zr, arb_midref(z), prec);
        else
            mag_swap(arb_radref(z), zr);

        mag_clear(ym);
        mag_clear(zr);

    }
}

void
arb_mul(arb_t z, const arb_t x, const arb_t y, long prec)
{
    mag_t zr, xm, ym;
    int inexact;

    if (arb_is_exact(x))
    {
        arb_mul_arf(z, y, arb_midref(x), prec);
    }
    else if (arb_is_exact(y))
    {
        arb_mul_arf(z, x, arb_midref(y), prec);
    }
    else if (ARB_IS_LAGOM(x) && ARB_IS_LAGOM(y) && ARB_IS_LAGOM(z))
    {
        mag_fast_init_set_arf(xm, arb_midref(x));
        mag_fast_init_set_arf(ym, arb_midref(y));

        mag_init(zr);
        mag_fast_mul(zr, xm, arb_radref(y));
        mag_fast_addmul(zr, ym, arb_radref(x));
        mag_fast_addmul(zr, arb_radref(x), arb_radref(y));

        inexact = arf_mul(arb_midref(z), arb_midref(x), arb_midref(y), prec, ARB_RND);

        if (inexact)
            arf_mag_fast_add_ulp(zr, zr, arb_midref(z), prec);

        *arb_radref(z) = *zr;
    }
    else
    {
        mag_init_set_arf(xm, arb_midref(x));
        mag_init_set_arf(ym, arb_midref(y));

        mag_init(zr);
        mag_mul(zr, xm, arb_radref(y));
        mag_addmul(zr, ym, arb_radref(x));
        mag_addmul(zr, arb_radref(x), arb_radref(y));

        inexact = arf_mul(arb_midref(z), arb_midref(x), arb_midref(y), prec, ARB_RND);

        if (inexact)
            arf_mag_add_ulp(arb_radref(z), zr, arb_midref(z), prec);
        else
            mag_swap(arb_radref(z), zr);

        mag_clear(xm);
        mag_clear(ym);
        mag_clear(zr);
    }
}

void
arb_mul_si(arb_t z, const arb_t x, long y, long prec)
{
    arf_t t;
    arf_init_set_si(t, y); /* no need to free */
    arb_mul_arf(z, x, t, prec);
}

void
arb_mul_ui(arb_t z, const arb_t x, ulong y, long prec)
{
    arf_t t;
    arf_init_set_ui(t, y); /* no need to free */
    arb_mul_arf(z, x, t, prec);
}

void
arb_mul_fmpz(arb_t z, const arb_t x, const fmpz_t y, long prec)
{
    arf_t t;

    if (!COEFF_IS_MPZ(*y))
    {
        arf_init_set_si(t, *y); /* no need to free */
        arb_mul_arf(z, x, t, prec);
    }
    else
    {
        arf_init(t);
        arf_set_fmpz(t, y);
        arb_mul_arf(z, x, t, prec);
        arf_clear(t);
    }
}

