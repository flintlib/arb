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

#include "fmprb.h"

void
fmprb_acos(fmprb_t z, const fmprb_t x, long prec)
{
    fmprb_t t, u;

    if (fmprb_is_exact(x))
    {
        int side;

        if (fmpr_is_zero(fmprb_midref(x)))
        {
            fmprb_const_pi(z, prec);
            fmprb_mul_2exp_si(z, z, -1);
            return;
        }

        side = fmpr_cmpabs_2exp_si(fmprb_midref(x), 0);

        /* +/- 1 */
        if (side == 0)
        {
            if (fmpr_is_one(fmprb_midref(x)))
                fmprb_zero(z);
            else
                fmprb_const_pi(z, prec);
            return;
        }
        else if (side > 0)
        {
            fmpr_nan(fmprb_midref(z));
            fmpr_pos_inf(fmprb_radref(z));
            return;
        }
    }

    fmprb_init(t);
    fmprb_init(u);

    fmprb_mul(t, x, x, FMPR_PREC_EXACT);
    fmprb_sub_ui(t, t, 1, prec);
    fmprb_neg(t, t);
    fmprb_rsqrt(t, t, prec);
    fmprb_mul(t, x, t, prec);
    fmprb_atan(t, t, prec);
    fmprb_const_pi(u, prec);
    fmprb_mul_2exp_si(u, u, -1);
    fmprb_sub(z, u, t, prec);

    fmprb_clear(t);
    fmprb_clear(u);
}

