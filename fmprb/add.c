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
fmprb_add(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec)
{
    long r;

    fmpr_add(fmprb_radref(z), fmprb_radref(x), fmprb_radref(y), FMPRB_RAD_PREC, FMPR_RND_UP);
    r = fmpr_add(fmprb_midref(z), fmprb_midref(x), fmprb_midref(y), prec, FMPR_RND_DOWN);

    fmpr_add_error_result(fmprb_radref(z), fmprb_radref(z), fmprb_midref(z),
        r, FMPRB_RAD_PREC, FMPR_RND_UP);

    fmprb_adjust(z);
}

void
fmprb_add_ui(fmprb_t z, const fmprb_t x, ulong y, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_set_ui(t, y);
    fmprb_add(z, x, t, prec);
    fmprb_clear(t);
}

void
fmprb_add_si(fmprb_t z, const fmprb_t x, long y, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_set_si(t, y);
    fmprb_add(z, x, t, prec);
    fmprb_clear(t);
}

void
fmprb_add_fmpz(fmprb_t z, const fmprb_t x, const fmpz_t y, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_set_fmpz(t, y);
    fmprb_add(z, x, t, prec);
    fmprb_clear(t);
}

void
fmprb_add_fmpr(fmprb_t z, const fmprb_t x, const fmpr_t y, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_set_fmpr(t, y);
    fmprb_add(z, x, t, prec);
    fmprb_clear(t);
}
