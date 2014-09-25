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
arb_get_interval_mpfr(mpfr_t a, mpfr_t b, const arb_t x)
{
    arf_t r, t;
    arf_init(t);
    arf_init_set_mag_shallow(r, arb_radref(x));
    arf_sub(t, arb_midref(x), r, mpfr_get_prec(a), ARF_RND_FLOOR);
    arf_get_mpfr(a, t, MPFR_RNDD);
    arf_add(t, arb_midref(x), r, mpfr_get_prec(b), ARF_RND_CEIL);
    arf_get_mpfr(b, t, MPFR_RNDU);
    arf_clear(t);
}

