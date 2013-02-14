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
gamma_taylor_bound_remainder(fmpr_t err, const fmpr_t z, long n)
{
    fmpr_t t, u;

    gamma_taylor_bound_extend_cache(n + 1);

    fmpr_init(t);
    fmpr_init(u);

    /* denominator: 1 - r(n+1) * z, rounded down */
    fmpr_mul(t, gamma_taylor_bound_ratio_cache + n + 1,
        z, FMPRB_RAD_PREC, FMPR_RND_UP);
    fmpr_one(u);
    fmpr_sub(u, u, t, FMPRB_RAD_PREC, FMPR_RND_DOWN);

    if (fmpr_sgn(u) <= 0)
    {
        fmpr_pos_inf(err);
    }
    else
    {
        fmpr_pow_sloppy_ui(t, z, n, FMPRB_RAD_PREC, FMPR_RND_UP);
        fmpr_mul_2exp_si(t, t, gamma_taylor_bound_mag_cache[n + 1]);
        fmpr_div(err, t, u, FMPRB_RAD_PREC, FMPR_RND_UP);
    }

    fmpr_clear(t);
    fmpr_clear(u);
}

