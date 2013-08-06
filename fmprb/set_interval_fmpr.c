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

#include "fmprb.h"

void
fmprb_set_interval_fmpr(fmprb_t x, const fmpr_t a, const fmpr_t b, long prec)
{
    fmpr_t t;
    fmpr_init(t);
    fmprb_set_fmpr(x, a);
    fmprb_add_fmpr(x, x, b, prec);
    fmpr_sub(t, b, a, FMPRB_RAD_PREC, FMPR_RND_UP);
    fmprb_add_error_fmpr(x, t);
    if (fmpr_is_nan(fmprb_radref(x)))
        fmpr_pos_inf(fmprb_radref(x));
    fmprb_mul_2exp_si(x, x, -1);
    fmpr_clear(t);
}

