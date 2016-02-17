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

    Copyright (C) 2016 Fredrik Johansson

******************************************************************************/

#include "arb.h"

void
arb_sgn(arb_t res, const arb_t x)
{
    if (arb_is_zero(x))
    {
        arb_zero(res);
    }
    else if (arb_contains_zero(x))
    {
        arf_zero(arb_midref(res));
        mag_one(arb_radref(res));
    }
    else
    {
        arb_set_si(res, arf_sgn(arb_midref(x)));
    }
}

