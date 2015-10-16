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

    Copyright (C) 2014 Fredrik Johansson

******************************************************************************/

#include "arb.h"

void
arb_rising_ui(arb_t y, const arb_t x, ulong n, long prec)
{
    if (n < FLINT_MAX(prec, 100))
    {
        arb_rising_ui_rec(y, x, n, prec);
    }
    else
    {
        arb_t t;
        arb_init(t);
        arb_add_ui(t, x, n, prec);
        arb_gamma(t, t, prec);
        arb_rgamma(y, x, prec);
        arb_mul(y, y, t, prec);
        arb_clear(t);
    }
}

void
arb_rising(arb_t y, const arb_t x, const arb_t n, long prec)
{
    if (arb_is_int(n) && arf_sgn(arb_midref(n)) >= 0 &&
        arf_cmpabs_ui(arb_midref(n), FLINT_MAX(prec, 100)) < 0)
    {
        arb_rising_ui_rec(y, x,
            arf_get_si(arb_midref(n), ARF_RND_DOWN), prec);
    }
    else
    {
        arb_t t;
        arb_init(t);
        arb_add(t, x, n, prec);
        arb_gamma(t, t, prec);
        arb_rgamma(y, x, prec);
        arb_mul(y, y, t, prec);
        arb_clear(t);
    }
}

