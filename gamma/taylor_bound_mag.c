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

long
gamma_taylor_bound_mag(long n)
{
    fmprb_t t, u;
    long v;

    fmprb_init(t);
    fmprb_init(u);

    /* (pi-1) n */
    fmprb_const_pi(t, FMPRB_RAD_PREC);
    fmprb_sub_ui(t, t, 1, FMPRB_RAD_PREC);
    fmprb_mul_ui(t, t, n, FMPRB_RAD_PREC);

    /* (3-5n) log(n/6) */
    fmprb_set_ui(u, n);
    fmprb_div_ui(u, u, 6, FMPRB_RAD_PREC);
    fmprb_log(u, u, FMPRB_RAD_PREC);
    fmprb_mul_si(u, u, 3 - 5*n, FMPRB_RAD_PREC);

    fmprb_add(t, t, u, FMPRB_RAD_PREC);

    /* divide by 6 log(2) */
    fmprb_log_ui(u, 2, FMPRB_RAD_PREC);
    fmprb_mul_ui(u, u, 6, FMPRB_RAD_PREC);
    fmprb_div(t, t, u, FMPRB_RAD_PREC);

    /* upper bound */
    fmpr_add(fmprb_midref(t), fmprb_midref(t), fmprb_radref(t),
        FMPRB_RAD_PREC, FMPR_RND_CEIL);

    v = fmpr_get_si(fmprb_midref(t), FMPR_RND_CEIL);

    fmprb_clear(t);
    fmprb_clear(u);

    return v;
}

