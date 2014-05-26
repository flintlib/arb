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

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include "arb.h"

void
arb_acos(arb_t z, const arb_t x, long prec)
{
    arb_t t, u;

    if (arb_is_exact(x))
    {
        int side;

        if (arf_is_zero(arb_midref(x)))
        {
            arb_const_pi(z, prec);
            arb_mul_2exp_si(z, z, -1);
            return;
        }

        side = arf_cmpabs_2exp_si(arb_midref(x), 0);

        /* +/- 1 */
        if (side == 0)
        {
            if (arf_is_one(arb_midref(x)))
                arb_zero(z);
            else
                arb_const_pi(z, prec);
            return;
        }
        else if (side > 0)
        {
            arb_indeterminate(z);
            return;
        }
    }

    arb_init(t);
    arb_init(u);

    arb_one(t);
    arb_submul(t, x, x, prec);
    arb_rsqrt(t, t, prec);
    arb_mul(t, x, t, prec);
    arb_atan(t, t, prec);
    arb_const_pi(u, prec);
    arb_mul_2exp_si(u, u, -1);
    arb_sub(z, u, t, prec);

    arb_clear(t);
    arb_clear(u);
}

