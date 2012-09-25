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

long
_fmpr_atan(fmpr_t y, const fmpr_t x, long prec, fmpr_rnd_t rnd)
{
    if (fmpr_is_zero(x) || fmpr_is_nan(x))
    {
        fmpr_set(y, x);
        return FMPR_RESULT_EXACT;
    }
    else
    {
        long r;
        CALL_MPFR_FUNC(r, mpfr_atan, y, x, prec, rnd);
        return r;
    }
}

void
fmprb_atan(fmprb_t z, const fmprb_t x, long prec)
{
    long r;

    if (fmprb_is_exact(x))
    {
        r = _fmpr_atan(fmprb_midref(z), fmprb_midref(x), prec, FMPR_RND_DOWN);
        fmpr_set_error_result(fmprb_radref(z), fmprb_midref(z), r);
    }
    else
    {
        fmpr_t t, u;

        fmpr_init(t);
        fmpr_init(u);

        if (fmpr_sgn(fmprb_midref(x)) >= 0)
        {
            fmpr_sub(t, fmprb_midref(x), fmprb_radref(x), FMPRB_RAD_PREC, FMPR_RND_DOWN);
        }
        else
        {
            fmpr_add(t, fmprb_midref(x), fmprb_radref(x), FMPRB_RAD_PREC, FMPR_RND_DOWN);
            fmpr_neg(t, t);
        }

        if (fmpr_sgn(t) > 0)
        {
            fmpr_mul(t, t, t, FMPRB_RAD_PREC, FMPR_RND_DOWN);
            fmpr_add_ui(t, t, 1UL, FMPRB_RAD_PREC, FMPR_RND_DOWN);
            fmpr_div(t, fmprb_radref(x), t, FMPRB_RAD_PREC, FMPR_RND_UP);
        }
        else
        {
            fmpr_set(t, fmprb_radref(x));
        }

        r = _fmpr_atan(fmprb_midref(z), fmprb_midref(x), prec, FMPR_RND_DOWN);
        fmpr_add_error_result(fmprb_radref(z), t, fmprb_midref(z), r,
            FMPRB_RAD_PREC, FMPR_RND_UP);

        fmpr_clear(t);
        fmpr_clear(u);
    }

    fmprb_adjust(z);
}
