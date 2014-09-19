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

void arb_div_2expm1_ui(arb_t z, const arb_t x, ulong n, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    arb_get_fmprb(t, x);
    fmprb_div_2expm1_ui(t, t, n, prec);
    arb_set_fmprb(z, t);
    fmprb_clear(t);
}

void arb_root(arb_t z, const arb_t x, ulong k, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    arb_get_fmprb(t, x);
    fmprb_root(t, t, k, prec);
    arb_set_fmprb(z, t);
    fmprb_clear(t);
}

