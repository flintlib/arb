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

static __inline__ void
_arb_get_mag_lower_nonnegative(mag_t z, const arf_t mid, const mag_t rad)
{
    if (mag_is_inf(rad) || arf_is_special(mid) || arf_sgn(mid) < 0)
    {
        mag_zero(z);
    }
    else if (mag_is_zero(rad))
    {
         arf_get_mag_lower(z, mid);
    }
    else
    {
        long shift, fix;

        shift = _fmpz_sub_small(MAG_EXPREF(mid), MAG_EXPREF(rad));

        /* mid < rad */
        if (shift < 0)
        {
            mag_zero(z);
        }
        else if (shift <= 1) /* can be cancellation */
        {
            arf_t t;
            arf_init(t);

            arf_set_mag(t, rad);

            arf_sub(t, mid, t, MAG_BITS, ARF_RND_DOWN);

            if (arf_sgn(t) <= 0)
                mag_zero(z);
            else
                arf_get_mag_lower(z, t);

            arf_clear(t);
        }
        else
        {
            mp_limb_t m;

            ARF_GET_TOP_LIMB(m, mid);
            m = m >> (FLINT_BITS - MAG_BITS);

            if (shift <= MAG_BITS)
                m = m - (MAG_MAN(rad) >> shift) - 1;
            else
                m = m - 1;

            fix = !(m >> (MAG_BITS - 1));
            m <<= fix;
            MAG_MAN(z) = m;
            _fmpz_add_fast(MAG_EXPREF(z), MAG_EXPREF(mid), -fix);
        }
    }
}

void
arb_get_mag_lower_nonnegative(mag_t z, const arb_t x)
{
    _arb_get_mag_lower_nonnegative(z, arb_midref(x), arb_radref(x));
}
