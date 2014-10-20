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

#include "acb_modular.h"

void
acb_modular_elliptic_p(acb_t r, const acb_t z, const acb_t tau, long prec)
{
    acb_struct t0[4], tz[4];
    acb_t t;
    int i;

    acb_init(t);

    for (i = 0; i < 4; i++)
    {
        acb_init(t0 + i);
        acb_init(tz + i);
    }

    acb_modular_theta(tz + 0, tz + 1, tz + 2, tz + 3, z, tau, prec);
    acb_zero(t);
    acb_modular_theta(t0 + 0, t0 + 1, t0 + 2, t0 + 3, t, tau, prec);

    acb_mul(t, t0 + 1, t0 + 2, prec);
    acb_mul(t, t, tz + 3, prec);
    acb_div(t, t, tz + 0, prec);
    acb_mul(t, t, t, prec);

    acb_pow_ui(t0 + 1, t0 + 1, 4, prec);
    acb_pow_ui(t0 + 2, t0 + 2, 4, prec);
    acb_add(r, t0 + 1, t0 + 2, prec);
    acb_div_ui(r, r, 3, prec);

    acb_sub(r, t, r, prec);

    acb_const_pi(t, prec);
    acb_mul(t, t, t, prec);
    acb_mul(r, r, t, prec);

    acb_clear(t);

    for (i = 0; i < 4; i++)
    {
        acb_clear(t0 + i);
        acb_clear(tz + i);
    }
}

