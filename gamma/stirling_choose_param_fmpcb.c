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

#include "gamma.h"

/* tuning factor */
#define GAMMA_STIRLING_BETA 0.23

void
gamma_stirling_choose_param_fmpcb(int * reflect, long * r, long * n,
    const fmpcb_t z, int use_reflect, long prec)
{
    if (fmpr_cmpabs_2exp_si(fmprb_midref(fmpcb_realref(z)), 40) > 0 ||
        fmpr_cmpabs_2exp_si(fmprb_midref(fmpcb_imagref(z)), 40) > 0)
    {
        if (use_reflect && fmpr_sgn(fmprb_midref(fmpcb_realref(z))) < 0)
            *reflect = 1;
        else
            *reflect = 0;

        *r = 0;
        *n = 1;
    }
    else
    {
        double x, y;

        x = fmpr_get_d(fmprb_midref(fmpcb_realref(z)), FMPR_RND_NEAR);
        y = fmpr_get_d(fmprb_midref(fmpcb_imagref(z)), FMPR_RND_NEAR);

        gamma_stirling_choose_param(reflect, r, n, x, y,
            GAMMA_STIRLING_BETA, use_reflect, prec);
    }
}

