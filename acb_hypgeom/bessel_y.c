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

#include "acb_hypgeom.h"

static void
phase(acb_t res, const arb_t re, const arb_t im)
{
    if (arb_is_nonnegative(re) || arb_is_negative(im))
    {
        acb_one(res);
    }
    else if (arb_is_negative(re) && arb_is_nonnegative(im))
    {
        acb_set_si(res, -3);
    }
    else
    {
        arb_zero(acb_imagref(res));
        /* -1 +/- 2 */
        arf_set_si(arb_midref(acb_realref(res)), -1);
        mag_one(arb_radref(acb_realref(res)));
        mag_mul_2exp_si(arb_radref(acb_realref(res)), arb_radref(acb_realref(res)), 1);
    }
}

void
acb_hypgeom_bessel_y(acb_t res, const acb_t nu, const acb_t z, long prec)
{
    acb_t t, u, v;

    acb_init(t);
    acb_init(u);
    acb_init(v);

    if (acb_is_int(nu))
    {
        int is_real = acb_is_real(nu) && acb_is_real(z)
            && arb_is_positive(acb_realref(z));

        acb_mul_onei(t, z);
        acb_hypgeom_bessel_k(t, nu, t, prec);
        acb_onei(u);
        acb_pow(u, u, nu, prec);
        acb_mul(t, t, u, prec);
        acb_const_pi(u, prec);
        acb_div(t, t, u, prec);
        acb_mul_2exp_si(t, t, 1);
        acb_neg(t, t);

        acb_hypgeom_bessel_j(u, nu, z, prec);
        phase(v, acb_realref(z), acb_imagref(z));
        acb_mul(u, u, v, prec);
        acb_mul_onei(u, u);

        acb_sub(res, t, u, prec);

        if (is_real)
            arb_zero(acb_imagref(res));
    }
    else
    {
        acb_sin_cos_pi(t, u, nu, prec);

        acb_hypgeom_bessel_j(v, nu, z, prec);
        acb_mul(v, v, u, prec);

        acb_neg(u, nu);
        acb_hypgeom_bessel_j(u, u, z, prec);
        acb_sub(v, v, u, prec);

        acb_div(res, v, t, prec);
    }

    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
}

