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

void
arb_mul(arb_t z, const arb_t x, const arb_t y, long prec)
{
    mag_t zr, xm, ym;
    int inexact;

    if (ARB_IS_LAGOM(x) && ARB_IS_LAGOM(y) && ARB_IS_LAGOM(z))
    {
        mag_fast_init_set_arf(xm, ARB_MIDREF(x));
        mag_fast_init_set_arf(ym, ARB_MIDREF(y));

        mag_init(zr);
        mag_fast_mul(zr, xm, ARB_RADREF(y));
        mag_fast_addmul(zr, ym, ARB_RADREF(x));
        mag_fast_addmul(zr, ARB_RADREF(x), ARB_RADREF(y));

        inexact = arf_mul(ARB_MIDREF(z), ARB_MIDREF(x), ARB_MIDREF(y),
            prec, ARF_RND_DOWN);

        if (inexact)
            mag_fast_add_2exp_si(zr, zr, ARF_EXP(ARB_MIDREF(z)) - prec);

        *ARB_RADREF(z) = *zr;
    }
    else
    {
        printf("big exponents and nans/infs not implemented yet!\n");
        abort();
    }
}

