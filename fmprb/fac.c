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

static void
fmprb_fac_ui_bsplit(fmprb_t x, ulong a, ulong b, long prec)
{
    if (b - a == 1)
    {
        fmprb_set_ui(x, b);
    }
    else
    {
        fmprb_t y;
        long m = a + (b - a) / 2;

        fmprb_init(y);
        fmprb_fac_ui_bsplit(x, a, m, prec);
        fmprb_fac_ui_bsplit(y, m, b, prec);
        fmprb_mul(x, x, y, prec);

        fmprb_clear(y);
    }
}

void
fmprb_fac_ui(fmprb_t x, ulong n, long prec)
{
    fmprb_fac_ui_bsplit(x, ULONG_MAX, FLINT_MAX(n, 1),
        prec + FLINT_BIT_COUNT(n));
}
