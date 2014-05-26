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

    Copyright (C) 2012-2014 Fredrik Johansson

******************************************************************************/

#include "arb.h"

void
arb_add(arb_t z, const arb_t x, const arb_t y, long prec)
{
    int inexact;

    inexact = arf_add(arb_midref(z), arb_midref(x), arb_midref(y), prec, ARB_RND);

    mag_add(arb_radref(z), arb_radref(x), arb_radref(y));
    if (inexact)
        arf_mag_add_ulp(arb_radref(z), arb_radref(z), arb_midref(z), prec);
}

void
arb_add_arf(arb_t z, const arb_t x, const arf_t y, long prec)
{
    int inexact;
    inexact = arf_add(arb_midref(z), arb_midref(x), y, prec, ARB_RND);
    if (inexact)
        arf_mag_add_ulp(arb_radref(z), arb_radref(x), arb_midref(z), prec);
    else
        mag_set(arb_radref(z), arb_radref(x));
}

void
arb_add_ui(arb_t z, const arb_t x, ulong y, long prec)
{
    int inexact;
    inexact = arf_add_ui(arb_midref(z), arb_midref(x), y, prec, ARB_RND);
    if (inexact)
        arf_mag_add_ulp(arb_radref(z), arb_radref(x), arb_midref(z), prec);
    else
        mag_set(arb_radref(z), arb_radref(x));
}

void
arb_add_si(arb_t z, const arb_t x, long y, long prec)
{
    int inexact;
    inexact = arf_add_si(arb_midref(z), arb_midref(x), y, prec, ARB_RND);
    if (inexact)
        arf_mag_add_ulp(arb_radref(z), arb_radref(x), arb_midref(z), prec);
    else
        mag_set(arb_radref(z), arb_radref(x));
}

void
arb_add_fmpz(arb_t z, const arb_t x, const fmpz_t y, long prec)
{
    int inexact;
    inexact = arf_add_fmpz(arb_midref(z), arb_midref(x), y, prec, ARB_RND);
    if (inexact)
        arf_mag_add_ulp(arb_radref(z), arb_radref(x), arb_midref(z), prec);
    else
        mag_set(arb_radref(z), arb_radref(x));
}

void
arb_add_fmpz_2exp(arb_t z, const arb_t x, const fmpz_t man, const fmpz_t exp, long prec)
{
    int inexact;
    inexact = arf_add_fmpz_2exp(arb_midref(z), arb_midref(x), man, exp, prec, ARB_RND);
    if (inexact)
        arf_mag_add_ulp(arb_radref(z), arb_radref(x), arb_midref(z), prec);
    else
        mag_set(arb_radref(z), arb_radref(x));
}

