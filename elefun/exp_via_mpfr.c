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

#include "elefun.h"

static __inline__ void
elefun_exp_fmpr_via_mpfr(fmprb_t z, const fmpr_t x, long prec)
{
    long r;
    r = fmpr_exp(fmprb_midref(z), x, prec, FMPR_RND_DOWN);
    fmpr_set_error_result(fmprb_radref(z), fmprb_midref(z), r);
}

void
elefun_exp_via_mpfr(fmprb_t z, const fmprb_t x, long prec)
{
    if (fmprb_is_exact(x))
    {
        elefun_exp_fmpr_via_mpfr(z, fmprb_midref(x), prec);
    }
    else
    {
        /* exp(a+b) - exp(a) = exp(a) * (exp(b)-1) <= b * exp(a+b) */
        fmpr_t t;
        fmpr_init(t);

        fmpr_add(t, fmprb_midref(x), fmprb_radref(x), FMPRB_RAD_PREC, FMPR_RND_CEIL);
        fmpr_exp(t, t, FMPRB_RAD_PREC, FMPR_RND_UP);
        fmpr_mul(t, t, fmprb_radref(x), FMPRB_RAD_PREC, FMPR_RND_UP);

        elefun_exp_fmpr_via_mpfr(z, fmprb_midref(x), prec);
        fmpr_add(fmprb_radref(z), fmprb_radref(z), t, FMPRB_RAD_PREC, FMPR_RND_UP);

        fmpr_clear(t);
    }
}

