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
fmprb_log(fmprb_t z, const fmprb_t x, long prec)
{
    long r;

    if (fmprb_contains_zero(x) || fmpr_sgn(fmprb_midref(x)) < 0)
    {
        fmpr_nan(fmprb_midref(z));
        fmpr_pos_inf(fmprb_radref(z));
        return;
    }

    if (fmprb_is_exact(x))
    {
        r = fmpr_log(fmprb_midref(z), fmprb_midref(x), prec, FMPR_RND_DOWN);
        fmpr_set_error_result(fmprb_radref(z), fmprb_midref(z), r);
    }
    else
    {
        /*
        Let the input be [a-b, a+b]. We require a > b >= 0 (otherwise the
        interval contains zero or a negative number and the logarithm is not
        defined). The error is largest at a-b, and we have

        log(a) - log(a-b) = log(1 + b/(a-b)).
        */
        fmpr_t err;
        fmpr_init(err);
        fmpr_sub(err, fmprb_midref(x), fmprb_radref(x), FMPRB_RAD_PREC, FMPR_RND_DOWN);
        fmpr_div(err, fmprb_radref(x), err, FMPRB_RAD_PREC, FMPR_RND_UP);
        fmpr_log1p(err, err, FMPRB_RAD_PREC, FMPR_RND_UP);

        r = fmpr_log(fmprb_midref(z), fmprb_midref(x), prec, FMPR_RND_DOWN);
        fmpr_add_error_result(fmprb_radref(z), err, fmprb_midref(z), r,
            FMPRB_RAD_PREC, FMPR_RND_UP);
        fmpr_clear(err);
    }

    fmprb_adjust(z);
}

void
fmprb_log_ui(fmprb_t z, ulong x, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_set_ui(t, x);
    fmprb_log(z, t, prec);
    fmprb_clear(t);
}

void
fmprb_log_fmpz(fmprb_t z, const fmpz_t x, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_set_fmpz(t, x);
    fmprb_log(z, t, prec);
    fmprb_clear(t);
}
