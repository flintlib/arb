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

#include "arb.h"

void
arb_bin_ui(arb_t x, const arb_t n, ulong k, long prec)
{
    if (k == 0)
    {
        arb_one(x);
    }
    else if (k == 1)
    {
        arb_set_round(x, n, prec);
    }
    else
    {
        arb_t t, u;

        arb_init(t);
        arb_init(u);

        arb_sub_ui(t, n, k - 1, prec);
        arb_rising_ui(t, t, k, prec);
        arb_fac_ui(u, k, prec);
        arb_div(x, t, u, prec);

        arb_clear(t);
        arb_clear(u);
    }
}

void
arb_bin_uiui(arb_t x, ulong n, ulong k, long prec)
{
    arb_t t;
    arb_init(t);
    arb_set_ui(t, n);
    arb_bin_ui(x, t, k, prec);
    arb_clear(t);
}

