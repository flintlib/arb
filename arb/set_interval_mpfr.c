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
arb_set_interval_mpfr(arb_t x, const mpfr_t a, const mpfr_t b, long prec)
{
    arf_t aa, bb;

    arf_init(aa);
    arf_init(bb);

    arf_set_mpfr(aa, a);
    arf_set_mpfr(bb, b);

    arb_set_interval_arf(x, aa, bb, prec);

    arf_clear(aa);
    arf_clear(bb);
}

