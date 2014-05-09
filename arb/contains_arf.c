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

int
arb_contains_arf(const arb_t x, const arf_t y)
{
    fmprb_t t;
    fmpr_t u;
    int result;

    fmprb_init(t);
    fmpr_init(u);

    arb_get_fmprb(t, x);
    arf_get_fmpr(u, y);

    result = fmprb_contains_fmpr(t, u);

    fmprb_clear(t);
    fmpr_clear(u);

    return result;
}

