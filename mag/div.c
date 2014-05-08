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

#include "mag.h"

void
mag_div(mag_t z, const mag_t x, const mag_t y)
{
    if (mag_is_special(x) || mag_is_special(y))
    {
        if (mag_is_zero(y) || mag_is_inf(x))
            mag_inf(z);
        else
            mag_zero(z);
    }
    else
    {
        mp_limb_t q;
        long fix;

#if FLINT_BITS == 64
        q = (MAG_MAN(x) << MAG_BITS) / MAG_MAN(y) + 1;
#else
        mp_limb_t hi, lo, r;
        lo = MAG_MAN(x) << MAG_BITS;
        hi = MAG_MAN(x) >> (FLINT_BITS - MAG_BITS);
        udiv_qrnnd(q, r, hi, lo, MAG_MAN(y));
        q += 1;
#endif
        fix = q >> MAG_BITS;
        q = (q >> fix) + fix;

        /* could overflow to next power of two */
        if (q >> MAG_BITS)
        {
            q >>= 1;
            fix += 1;
        }

        MAG_MAN(z) = q;
        _fmpz_sub2_fast(MAG_EXPREF(z), MAG_EXPREF(x), MAG_EXPREF(y), fix);
    }
}

