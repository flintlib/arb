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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "fmprb.h"

void
fmprb_rsqrt(fmprb_t z, const fmprb_t x, long prec)
{
    long r;

    if (fmprb_contains_nonpositive(x))
    {
        fmpr_nan(fmprb_midref(z));
        fmpr_pos_inf(fmprb_radref(z));
        return;
    }

    if (fmprb_is_exact(x))
    {
        r = fmpr_rsqrt(fmprb_midref(z), fmprb_midref(x), prec, FMPR_RND_DOWN);
        fmpr_set_error_result(fmprb_radref(z), fmprb_midref(z), r);
    }
    else
    {
        fmpr_t t, err;

        fmpr_init(t);
        fmpr_init(err);

        /* lower point */
        fmpr_sub(t, fmprb_midref(x), fmprb_radref(x), FMPRB_RAD_PREC, FMPR_RND_DOWN);

        /* error bound: (1/2) t^(-3/2) * rad */
        fmpr_rsqrt(err, t, FMPRB_RAD_PREC, FMPR_RND_UP);
        fmpr_divappr_abs_ubound(err, err, t, FMPRB_RAD_PREC);
        fmpr_mul(err, err, fmprb_radref(x), FMPRB_RAD_PREC, FMPR_RND_UP);
        fmpr_mul_2exp_si(err, err, -1);

        r = fmpr_rsqrt(fmprb_midref(z), fmprb_midref(x), prec, FMPR_RND_DOWN);
        fmpr_add_error_result(fmprb_radref(z), err, fmprb_midref(z), r,
            FMPRB_RAD_PREC, FMPR_RND_UP);

        fmpr_clear(t);
        fmpr_clear(err);
    }

    fmprb_adjust(z);
}

void
fmprb_rsqrt_ui(fmprb_t z, ulong x, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_set_ui(t, x);
    fmprb_rsqrt(z, t, prec);
    fmprb_clear(t);
}

