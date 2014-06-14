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

#include "arb.h"

void
arb_zeta_ui_asymp(arb_t x, ulong s, long prec)
{
    arb_one(x);

    if (s != 2 && s > prec)
    {
        arb_add_error_2exp_si(x, -prec);
    }
    else
    {
        arf_t t;
        arf_init(t);
        arf_set_ui_2exp_si(t, 1, -s);
        arb_add_arf(x, x, t, prec);
        arb_add_error_2exp_si(x, 2 - (3 * s) / 2);
        arf_clear(t);
    }
}
