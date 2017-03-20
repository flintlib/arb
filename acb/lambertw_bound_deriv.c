/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
acb_lambertw_bound_deriv(mag_t res, const acb_t z, const acb_t ez1, const fmpz_t k)
{
    mag_t t, u;

    mag_init(t);
    mag_init(u);

    /* Main approximation: W'(z) ~ 1/|z|  */
    if (fmpz_is_zero(k))
    {
        /* 1/(1+|z|) */
        acb_get_mag_lower(res, z);
        mag_one(t);
        mag_add_lower(res, res, t);
        mag_div(res, t, res);
    }
    else if (fmpz_is_pm1(k))
    {
        /* 2/|z| */
        mag_set_ui(t, 2);
        acb_get_mag_lower(res, z);
        mag_div(res, t, res);
    }
    else
    {
        /* |W'(z)|  = |1/z| |W(z)/(1+W(z))|
                   <= |1/z| (2pi) / (2pi-1) */
        mag_set_ui_2exp_si(t, 77, -6);
        acb_get_mag_lower(res, z);
        mag_div(res, t, res);
    }

    /* Compute correction near the branch point */
    if (fmpz_is_zero(k) || fmpz_equal_si(k,-1) ||
        (fmpz_is_one(k) && !arb_is_nonnegative(acb_imagref(z))))
    {
        acb_t b;
        acb_init(b);

        if (fmpz_is_zero(k))
        {
            /* [-4,1] + [-2,2]i */
            arf_set_si_2exp_si(arb_midref(acb_realref(b)), -3, -1);
            mag_set_ui_2exp_si(arb_radref(acb_realref(b)), 5, -1);
            arf_zero(arb_midref(acb_imagref(b)));
            mag_set_ui(arb_radref(acb_imagref(b)), 2);
        }
        else
        {
            /* k = 1   [-1/2,-1/4] + [-1/8,0) i */
            /* k = -1  [-1/2,-1/4] + [0,1/8] i */
            arf_set_si_2exp_si(arb_midref(acb_realref(b)), -3, -3);
            mag_set_ui_2exp_si(arb_radref(acb_realref(b)), 1, -3);
            if (fmpz_is_one(k))
                arf_set_si_2exp_si(arb_midref(acb_imagref(b)), -1, -4);
            else
                arf_set_ui_2exp_si(arb_midref(acb_imagref(b)), 1, -4);
            mag_set_ui_2exp_si(arb_radref(acb_imagref(b)), 1, -4);
        }

        if (acb_overlaps(z, b))
        {
            /* add 2/sqrt(|ez+1|) */
            acb_get_mag_lower(u, ez1);
            mag_rsqrt(u, u);
            mag_mul_2exp_si(u, u, 1);
            mag_add(res, res, u);
        }

        acb_clear(b);
    }

    mag_clear(t);
    mag_clear(u);
}

