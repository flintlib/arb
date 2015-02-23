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

void
acb_hypgeom_erfc(acb_t res, const acb_t z, long prec)
{
    if (arb_is_positive(acb_realref(z)))
    {
        acb_t t, u;

        acb_init(t);
        acb_init(u);

        acb_set_si(t, 1);
        acb_mul_2exp_si(t, t, -1);
        acb_mul(u, z, z, prec);

        acb_hypgeom_gamma_upper(res, t, u, 0, prec);

        arb_const_sqrt_pi(acb_realref(t), prec);
        acb_div_arb(res, res, acb_realref(t), prec);

        acb_clear(t);
        acb_clear(u);
    }
    else
    {
        acb_hypgeom_erf(res, z, prec);
        acb_sub_ui(res, res, 1, prec);
        acb_neg(res, res);
    }
}

