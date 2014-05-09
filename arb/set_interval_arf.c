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
arb_set_interval_arf(arb_t x, const arf_t a, const arf_t b, long prec)
{
    fmprb_t t;
    fmpr_t u, v;

    fmprb_init(t);
    fmpr_init(u);
    fmpr_init(v);

    arf_get_fmpr(u, a);
    arf_get_fmpr(v, b);

    fmprb_set_interval_fmpr(t, u, v, prec);

    arb_set_fmprb(x, t);

    fmprb_clear(t);
    fmpr_clear(u);
    fmpr_clear(v);
}

