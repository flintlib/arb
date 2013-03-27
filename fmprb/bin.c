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
fmprb_bin_ui(fmprb_t x, const fmprb_t n, ulong k, long prec)
{
    if (k == 0)
    {
        fmprb_one(x);
    }
    else if (k == 1)
    {
        fmprb_set_round(x, n, prec);
    }
    else
    {
        fmprb_t t, u;

        fmprb_init(t);
        fmprb_init(u);

        fmprb_sub_ui(t, n, k - 1, prec);
        fmprb_rising_ui(t, t, k, prec);
        fmprb_fac_ui(u, k, prec);
        fmprb_div(x, t, u, prec);

        fmprb_clear(t);
        fmprb_clear(u);
    }
}

void
fmprb_bin_uiui(fmprb_t x, ulong n, ulong k, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_set_ui(t, n);
    fmprb_bin_ui(x, t, k, prec);
    fmprb_clear(t);
}

