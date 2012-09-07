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
fmprb_exp(fmprb_t z, const fmprb_t x, long prec)
{
    long r;

    if (fmprb_is_exact(x))
    {
        r = fmpr_exp(fmprb_midref(z), fmprb_midref(x), prec, FMPR_RND_DOWN);
        fmpr_set_error_result(fmprb_radref(z), fmprb_midref(z), r);
    }
    else
    {
        /*
        exp(a+b) - exp(a) = exp(a) * (exp(b)-1)
        */
        fmpr_t t, u;

        fmpr_init(t);
        fmpr_init(u);

        fmpr_exp(t, fmprb_midref(x), FMPRB_RAD_PREC, FMPR_RND_UP);
        fmpr_expm1(u, fmprb_radref(x), FMPRB_RAD_PREC, FMPR_RND_UP);
        fmpr_mul(t, t, u, FMPRB_RAD_PREC, FMPR_RND_UP);

        r = fmpr_exp(fmprb_midref(z), fmprb_midref(x), prec, FMPR_RND_DOWN);
        fmpr_add_error_result(fmprb_radref(z), t, fmprb_midref(z), r,
            FMPRB_RAD_PREC, FMPR_RND_UP);

        fmpr_clear(t);
        fmpr_clear(u);
    }

    fmprb_adjust(z);
}
