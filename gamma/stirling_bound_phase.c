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

void
gamma_stirling_bound_phase(fmpr_t bound, const fmpcb_t z, long prec)
{
    fmpr_t x, y, t, u;
    int xsign;

    fmpr_init(x);
    fmpr_init(y);
    fmpr_init(t);
    fmpr_init(u);

    /* first compute x, y such that |arg(z)| <= arg(x+yi) */

    /* argument increases with smaller real parts */
    fmpr_sub(x, fmprb_midref(fmpcb_realref(z)),
        fmprb_radref(fmpcb_realref(z)), prec, FMPR_RND_FLOOR);

    xsign = fmpr_sgn(x);

    if (xsign >= 0)  /* argument increases away from the real axis */
        fmprb_get_abs_ubound_fmpr(y, fmpcb_imagref(z), prec);
    else  /* argument increases closer to the real axis */
        fmprb_get_abs_lbound_fmpr(y, fmpcb_imagref(z), prec);

    if (fmpr_is_zero(y))
    {
        if (xsign > 0)
            fmpr_one(bound);
        else
            fmpr_pos_inf(bound);
    }
    else
    {
        if (xsign >= 0)
        {
            /* compute upper bound for t = y / (sqrt(x^2 + y^2) + x) */
            fmpr_mul(t, x, x, prec, FMPR_RND_DOWN);
            fmpr_mul(u, y, y, prec, FMPR_RND_DOWN);
            fmpr_add(t, t, u, prec, FMPR_RND_DOWN);
            fmpr_sqrt(t, t, prec, FMPR_RND_DOWN);
            fmpr_add(t, t, x, prec, FMPR_RND_DOWN);
            fmpr_div(t, y, t, prec, FMPR_RND_UP);
        }
        else
        {
            /* compute upper bound for t = (sqrt(x^2 + y^2) - x) / y */
            fmpr_mul(t, x, x, prec, FMPR_RND_UP);
            fmpr_mul(u, y, y, prec, FMPR_RND_UP);
            fmpr_add(t, t, u, prec, FMPR_RND_UP);
            fmpr_sqrt(t, t, prec, FMPR_RND_UP);
            fmpr_sub(t, t, x, prec, FMPR_RND_UP);
            fmpr_div(t, t, y, prec, FMPR_RND_UP);
        }

        /* compute upper bound for sqrt(1 + t^2) */
        fmpr_mul(t, t, t, prec, FMPR_RND_UP);
        fmpr_add_ui(t, t, 1, prec, FMPR_RND_UP);
        fmpr_sqrt(bound, t, prec, FMPR_RND_UP);
    }

    fmpr_clear(x);
    fmpr_clear(y);
    fmpr_clear(t);
    fmpr_clear(u);
}

