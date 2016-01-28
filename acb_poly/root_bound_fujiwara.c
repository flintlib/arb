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

#include "acb_poly.h"

void
_acb_poly_root_bound_fujiwara(mag_t bound, acb_srcptr poly, slong len)
{
    mag_t t, u, v;
    slong i;

    if (len <= 1)
    {
        mag_inf(bound);
        return;
    }

    mag_init(t);
    mag_init(u);
    mag_init(v);

    /* u = 1/leading */
    acb_get_mag_lower(t, poly + len - 1);
    mag_one(u);
    mag_div(u, u, t);

    mag_zero(v);

    for (i = 0; i < len - 1; i++)
    {
        acb_get_mag(t, poly + len - 2 - i);
        mag_mul(t, t, u);

        if (i == len - 2)
            mag_mul_2exp_si(t, t, -1);

        mag_root(t, t, i + 1);
        mag_max(v, v, t);
    }

    mag_mul_2exp_si(bound, v, 1);

    mag_clear(t);
    mag_clear(u);
    mag_clear(v);
}

void
acb_poly_root_bound_fujiwara(mag_t bound, acb_poly_t poly)
{
    _acb_poly_root_bound_fujiwara(bound, poly->coeffs, poly->length);
}

