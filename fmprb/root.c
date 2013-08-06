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

    Copyright (C) 2012, 2013 Fredrik Johansson

******************************************************************************/

#include "fmprb.h"

void
fmprb_root(fmprb_t z, const fmprb_t x, ulong k, long prec)
{
    long r;

    if (k == 1)
    {
        fmprb_set_round(z, x, prec);
        return;
    }

    if (k == 2)
    {
        fmprb_sqrt(z, x, prec);
        return;
    }

    if (k == -2)
    {
        fmprb_rsqrt(z, x, prec);
        return;
    }

    if (fmprb_contains_negative(x))
    {
        fmprb_indeterminate(z);
        return;
    }

    if (fmprb_is_exact(x))
    {
        r = fmpr_root(fmprb_midref(z), fmprb_midref(x), k, prec, FMPR_RND_DOWN);
        fmpr_set_error_result(fmprb_radref(z), fmprb_midref(z), r);
    }
    else
    {
        fmpr_t t, u, err;

        fmpr_init(t);
        fmpr_init(u);
        fmpr_init(err);

        /* lower point */
        fmpr_sub(err, fmprb_midref(x), fmprb_radref(x), FMPRB_RAD_PREC, FMPR_RND_DOWN);

        if (fmpr_is_zero(err))
        {
            fmpr_root(err, fmprb_midref(x), k, FMPRB_RAD_PREC, FMPR_RND_UP);
            r = fmpr_root(fmprb_midref(z), fmprb_midref(x), k, prec, FMPR_RND_DOWN);
            fmpr_add_error_result(fmprb_radref(z), err, fmprb_midref(z), r,
                FMPRB_RAD_PREC, FMPR_RND_UP);
        }
        else
        {
            /* derivative x^(1/k) / (x k) at lower point */
            fmpr_root(t, err, k, FMPRB_RAD_PREC, FMPR_RND_UP);
            fmpr_mul_ui(u, err, k, FMPRB_RAD_PREC, FMPR_RND_DOWN);
            fmpr_divappr_abs_ubound(t, t, u, FMPRB_RAD_PREC);

            /* multiplied by distance */
            fmpr_mul(err, t, fmprb_radref(x), FMPRB_RAD_PREC, FMPR_RND_UP);

            r = fmpr_root(fmprb_midref(z), fmprb_midref(x), k, prec, FMPR_RND_DOWN);
            fmpr_add_error_result(fmprb_radref(z), err, fmprb_midref(z), r,
                FMPRB_RAD_PREC, FMPR_RND_UP);
        }

        fmpr_clear(t);
        fmpr_clear(u);
        fmpr_clear(err);
    }

    fmprb_adjust(z);
}

