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

void mag_add_arf_ulp(mag_t r, const arf_t z, long prec)
{
    if (ARF_IS_SPECIAL(z))
    {
        printf("error: ulp error not defined for special value!\n");
        abort();
    }
    else
    {
        /* todo: speed up case r is special here */
        fmpz_t e;
        fmpz_init(e);
        _fmpz_add_fast(e, ARF_EXPREF(z), -prec);
        mag_add_2exp_fmpz(r, r, e);
        fmpz_clear(e);
    }
}

void
mag_print(const mag_t x)
{
    fmpr_t t;
    fmpr_init(t);
    mag_get_fmpr(t, x);
    fmpr_printd(t, 5);
    fmpr_clear(t);
}

void
arb_addmul_slow(arb_t z, const arb_t x, const arb_t y, long prec)
{
    mag_t zr, xm, ym;
    int inexact;

    mag_init_set_arf(xm, ARB_MIDREF(x));
    mag_init_set_arf(ym, ARB_MIDREF(y));

    mag_init_set(zr, ARB_RADREF(z));
    mag_addmul(zr, xm, ARB_RADREF(y));
    mag_addmul(zr, ym, ARB_RADREF(x));
    mag_addmul(zr, ARB_RADREF(x), ARB_RADREF(y));

    inexact = arf_addmul(ARB_MIDREF(z), ARB_MIDREF(x), ARB_MIDREF(y),
        prec, ARF_RND_DOWN);

    if (inexact)
        mag_add_arf_ulp(zr, ARB_MIDREF(z), prec);

    mag_set(ARB_RADREF(z), zr); /* swap? */

    mag_clear(zr);
    mag_clear(xm);
    mag_clear(ym);
}

void
arb_addmul(arb_t z, const arb_t x, const arb_t y, long prec)
{
    mag_t zr, xm, ym;
    int inexact;

    if (ARB_IS_LAGOM(x) && ARB_IS_LAGOM(y) && ARB_IS_LAGOM(z))
    {
        mag_fast_init_set_arf(xm, ARB_MIDREF(x));
        mag_fast_init_set_arf(ym, ARB_MIDREF(y));

        mag_fast_init_set(zr, ARB_RADREF(z));
        mag_fast_addmul(zr, xm, ARB_RADREF(y));
        mag_fast_addmul(zr, ym, ARB_RADREF(x));
        mag_fast_addmul(zr, ARB_RADREF(x), ARB_RADREF(y));

        inexact = arf_addmul(ARB_MIDREF(z), ARB_MIDREF(x), ARB_MIDREF(y),
            prec, ARF_RND_DOWN);

        if (inexact)
            mag_fast_add_2exp_si(zr, zr, ARF_EXP(ARB_MIDREF(z)) - prec);

        *ARB_RADREF(z) = *zr;
    }
    else
    {
        arb_addmul_slow(z, x, y, prec);
    }
}

