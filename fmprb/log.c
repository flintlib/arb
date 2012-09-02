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

    if (!fmprb_is_exact(x))
        abort();

    r = fmpr_log(fmprb_midref(z), fmprb_midref(x), prec, FMPR_RND_DOWN);

    fmpr_set_error_result(fmprb_radref(z), fmprb_midref(z), r);

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
