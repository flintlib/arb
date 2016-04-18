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

    Copyright (C) 2016 Arb authors

******************************************************************************/

#include "acb_hypgeom.h"

void
acb_hypgeom_gamma_lower(acb_t res, const acb_t s, const acb_t z, int regularized, slong prec)
{
    acb_t s1, nz, t;

    if (!acb_is_finite(s) || !acb_is_finite(z))
    {
        acb_indeterminate(res);
        return;
    }

    acb_init(s1);
    acb_init(nz);
    acb_init(t);

    acb_add_ui(s1, s, 1, prec);
    acb_neg(nz, z);

    if (regularized == 0)
    {
        /* \gamma(s, z) = s^-1 z^s 1F1(s, 1+s, -z) */
        acb_hypgeom_m(res, s, s1, nz, 0, prec);
        acb_pow(t, z, s, prec);
        acb_mul(res, res, t, prec);
        acb_div(res, res, s, prec);
    }
    else if (regularized == 1)
    {
        /* P(s, z) = z^s \gamma^{*}(s, z) */
        acb_hypgeom_m(res, s, s1, nz, 1, prec);
        acb_pow(t, z, s, prec);
        acb_mul(res, res, t, prec);
    }
    else if (regularized == 2)
    {
        /* \gamma^{*}(s, z) */
        acb_hypgeom_m(res, s, s1, nz, 1, prec);
    }

    acb_clear(s1);
    acb_clear(nz);
    acb_clear(t);
}
