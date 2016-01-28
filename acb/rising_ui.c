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

#include "acb.h"

void
acb_rising_ui(acb_t y, const acb_t x, ulong n, slong prec)
{
    if (n < FLINT_MAX(prec, 100))
    {
        acb_rising_ui_rec(y, x, n, prec);
    }
    else
    {
        acb_t t;
        acb_init(t);
        acb_add_ui(t, x, n, prec);
        acb_gamma(t, t, prec);
        acb_rgamma(y, x, prec);
        acb_mul(y, y, t, prec);
        acb_clear(t);
    }
}

void
acb_rising(acb_t y, const acb_t x, const acb_t n, slong prec)
{
    if (acb_is_int(n) && arf_sgn(arb_midref(acb_realref(n))) >= 0 &&
        arf_cmpabs_ui(arb_midref(acb_realref(n)), FLINT_MAX(prec, 100)) < 0)
    {
        acb_rising_ui_rec(y, x,
            arf_get_si(arb_midref(acb_realref(n)), ARF_RND_DOWN), prec);
    }
    else
    {
        acb_t t;
        acb_init(t);
        acb_add(t, x, n, prec);
        acb_gamma(t, t, prec);
        acb_rgamma(y, x, prec);
        acb_mul(y, y, t, prec);
        acb_clear(t);
    }
}

