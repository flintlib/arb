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
#include "fmprb_poly.h"

void
gamma_taylor_fmprb(fmprb_t y, const fmprb_t x, long prec)
{
    long v;
    fmprb_t t, u;

    /* nearest integer (TODO: clamp, avoiding overflow) */
    v = fmpr_get_si(fmprb_midref(x), FMPR_RND_NEAR);

    fmprb_init(t);
    fmprb_init(u);

    if (v == 0)
    {
        gamma_taylor_eval_fmprb(u, x, prec);
        fmprb_mul(u, u, x, prec);
        fmprb_inv(y, u, prec);
    }
    else if (v > 0)
    {
        fmprb_sub_si(t, x, v, prec);
        gamma_taylor_eval_fmprb(t, t, prec);
        fmprb_sub_si(u, x, v - 1, prec);
        gamma_rising_fmprb_ui_bsplit(u, u, v - 1, prec);
        fmprb_div(y, u, t, prec);
    }
    else
    {
        fmprb_add_si(t, x, (-v), prec);
        gamma_taylor_eval_fmprb(t, t, prec);
        gamma_rising_fmprb_ui_bsplit(u, x, (-v) + 1, prec);
        fmprb_mul(y, u, t, prec);
        fmprb_inv(y, y, prec);
    }

    fmprb_clear(t);
    fmprb_clear(u);
}

