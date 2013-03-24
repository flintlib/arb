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
fmprb_add_error_fmpr(fmprb_t x, const fmpr_t err)
{
    fmpr_add(fmprb_radref(x), fmprb_radref(x), err, FMPRB_RAD_PREC, FMPR_RND_UP);
}

void
fmprb_add_error_2exp_si(fmprb_t x, long err)
{
    fmpr_t t;
    fmpr_init(t);
    fmpr_set_ui_2exp_si(t, 1, err);
    fmprb_add_error_fmpr(x, t);
    fmpr_clear(t);
}

void
fmprb_add_error_2exp_fmpz(fmprb_t x, const fmpz_t err)
{
    fmpr_t t;
    fmpr_init(t);
    fmpz_one(fmpr_manref(t));
    fmpz_set(fmpr_expref(t), err);
    fmprb_add_error_fmpr(x, t);
    fmpr_clear(t);
}

void
fmprb_add_error(fmprb_t x, const fmprb_t error)
{
    fmpr_t high;
    fmpr_init(high);

    fmpr_add(high, fmprb_midref(error), fmprb_radref(error),
        FMPRB_RAD_PREC, FMPR_RND_UP);

    fmpr_add(fmprb_radref(x), fmprb_radref(x), high,
        FMPRB_RAD_PREC, FMPR_RND_UP);

    fmpr_clear(high);
}
