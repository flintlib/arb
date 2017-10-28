/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_sinh(arb_t s, const arb_t x, slong prec)
{
    if (arb_is_zero(x))
    {
        arb_zero(s);
    }
    else if (!arb_is_finite(x))
    {
        if (arf_is_nan(arb_midref(x)))
            arb_indeterminate(s);
        else
            arb_zero_pm_inf(s);
    }
    else
    {
        arb_t t;
        slong wp = prec + 4;

        arb_init(t);

        if (arf_cmpabs_2exp_si(arb_midref(x), -1) <= 0 &&
               mag_cmp_2exp_si(arb_radref(x), -4) <= 0)
        {
            arb_expm1(s, x, wp);
            arb_add_ui(t, s, 1, wp);
            arb_div(t, s, t, wp);
            arb_add(s, s, t, prec);
        }
        else
        {
            arb_exp_invexp(s, t, x, wp);
            arb_sub(s, s, t, prec);
        }

        arb_mul_2exp_si(s, s, -1);
        arb_clear(t);
    }
}

void
arb_cosh(arb_t c, const arb_t x, slong prec)
{
    if (arb_is_zero(x))
    {
        arb_one(c);
    }
    else if (!arb_is_finite(x))
    {
        if (arf_is_nan(arb_midref(x)))
            arb_indeterminate(c);
        else
            arb_zero_pm_inf(c);
    }
    else
    {
        arb_t t;
        slong wp = prec + 4;

        arb_init(t);

        arb_exp_invexp(c, t, x, wp);
        arb_add(c, c, t, prec);
        arb_mul_2exp_si(c, c, -1);

        arb_clear(t);
    }
}

void
arb_sinh_cosh(arb_t s, arb_t c, const arb_t x, slong prec)
{
    if (arb_is_zero(x))
    {
        arb_zero(s);
        arb_one(c);
    }
    else if (!arb_is_finite(x))
    {
        if (arf_is_nan(arb_midref(x)))
        {
            arb_indeterminate(s);
            arb_indeterminate(c);
        }
        else
        {
            arb_zero_pm_inf(s);
            arb_zero_pm_inf(c);
        }
    }
    else
    {
        slong wp = prec + 4;

        arb_t t;
        arb_init(t);

        if (arf_cmpabs_2exp_si(arb_midref(x), -1) <= 0 &&
            mag_cmp_2exp_si(arb_radref(x), -4) <= 0)
        {
            arb_expm1(s, x, wp);
            arb_add_ui(t, s, 1, wp);
            arb_inv(c, t, wp);
            arb_addmul(s, s, c, prec);
            arb_add(c, c, t, prec);
        }
        else
        {
            arb_exp_invexp(c, t, x, wp);
            arb_sub(s, c, t, prec);
            arb_add(c, c, t, prec);
        }

        arb_mul_2exp_si(s, s, -1);
        arb_mul_2exp_si(c, c, -1);
        arb_clear(t);
    }
}

