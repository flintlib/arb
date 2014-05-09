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
arb_union(arb_t z, const arb_t x, const arb_t y, long prec)
{
    fmprb_t t, u, v;

    fmprb_init(t);
    fmprb_init(u);
    fmprb_init(v);

    arb_get_fmprb(t, x);
    arb_get_fmprb(u, y);

    fmprb_union(v, t, u, prec);

    arb_set_fmprb(z, v);

    fmprb_clear(t);
    fmprb_clear(u);
    fmprb_clear(v);
}

