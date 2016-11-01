/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_div_arf(arb_t z, const arb_t x, const arf_t y, slong prec)
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
arb_div(arb_t z, const arb_t x, const arb_t y, slong prec)
{
    int inexact;

    if (arb_is_exact(y))
    {
        arb_div_arf(z, x, arb_midref(y), prec);
    }
    else if (mag_is_inf(arb_radref(x)) || mag_is_inf(arb_radref(y)) || arf_is_zero(arb_midref(y)))
    {
        arf_div(arb_midref(z), arb_midref(x), arb_midref(y), prec, ARB_RND);
        mag_inf(arb_radref(z));
    }
    else if (ARB_IS_LAGOM(x) && ARB_IS_LAGOM(y) &&  /* fast path */
            MAG_EXP(arb_radref(y)) < ARF_EXP(arb_midref(y)) - 20 && prec > 20)
    {
        mag_t t, u, v;
        mag_init(t);  /* no need to free */
        mag_init(u);
        mag_init(v);

        /*               (x*yrad + y*xrad)/(y*(y-yrad))
             <=  (1+eps) (x*yrad + y*xrad)/y^2
             <=  (1+eps) (x*yrad/y^2 + y*xrad/y^2)
             <=  (1+eps) ((x/y)*yrad/y + xrad/y)
             <=  (1+eps) ((x/y)*yrad + xrad)/y
        */

        /* t = y */
        arf_get_mag_lower(t, arb_midref(y));

        inexact = arf_div(arb_midref(z), arb_midref(x), arb_midref(y), prec, ARB_RND);

        /* u = x/y */
        arf_get_mag(u, arb_midref(z));

        /* v = (x/y)*yrad + xrad */
        MAG_MAN(v) = MAG_MAN(arb_radref(x));
        MAG_EXP(v) = MAG_EXP(arb_radref(x));
        mag_fast_addmul(v, arb_radref(y), u);
        /* v = v / y */
        mag_div(arb_radref(z), v, t);

        /* correct for errors */
        MAG_MAN(t) = MAG_ONE_HALF + (MAG_ONE_HALF >> 16);
        MAG_EXP(t) = 1;
        mag_fast_mul(arb_radref(z), arb_radref(z), t);

        if (inexact)
            arf_mag_fast_add_ulp(arb_radref(z), arb_radref(z), arb_midref(z), prec);
    }
    else
    {
        mag_t zr, xm, ym, yl, yw;

        mag_init_set_arf(xm, arb_midref(x));
        mag_init_set_arf(ym, arb_midref(y));
        mag_init(zr);
        mag_init(yl);
        mag_init(yw);

        /* (|x|*yrad + |y|*xrad)/(|y|*(|y|-yrad)) */
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
arb_div_si(arb_t z, const arb_t x, slong y, slong prec)
{
    arf_t t;
    arf_init_set_si(t, y); /* no need to free */
    arb_div_arf(z, x, t, prec);
}

void
arb_div_ui(arb_t z, const arb_t x, ulong y, slong prec)
{
    arf_t t;
    arf_init_set_ui(t, y); /* no need to free */
    arb_div_arf(z, x, t, prec);
}

void
arb_div_fmpz(arb_t z, const arb_t x, const fmpz_t y, slong prec)
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
arb_fmpz_div_fmpz(arb_t z, const fmpz_t x, const fmpz_t y, slong prec)
{
    int inexact;

    inexact = arf_fmpz_div_fmpz(arb_midref(z), x, y, prec, ARB_RND);

    if (inexact)
        arf_mag_set_ulp(arb_radref(z), arb_midref(z), prec);
    else
        mag_zero(arb_radref(z));
}

void
arb_ui_div(arb_t z, ulong x, const arb_t y, slong prec)
{
    arb_t t;
    arb_init(t);
    arb_set_ui(t, x);
    arb_div(z, t, y, prec);
    arb_clear(t);
}

