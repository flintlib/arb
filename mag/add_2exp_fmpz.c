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
mag_add_2exp_fmpz(mag_t z, const mag_t x, const fmpz_t e)
{
    if (mag_is_special(x))
    {
        if (mag_is_zero(x))
        {
            MAG_MAN(z) = MAG_ONE_HALF;
            _fmpz_add_fast(MAG_EXPREF(z), e, 1);
        }
        else
        {
            mag_inf(z);
        }
    }
    else
    {
        long shift;

        shift = _fmpz_sub_small(MAG_EXPREF(x), e);

        if (shift > 0)
        {
            _fmpz_set_fast(MAG_EXPREF(z), MAG_EXPREF(x));

            if (shift >= MAG_BITS)
                MAG_MAN(z) = MAG_MAN(x) + LIMB_ONE;
            else
                MAG_MAN(z) = MAG_MAN(x) + (LIMB_ONE << (MAG_BITS - shift));
        }
        else
        {
            shift = -shift;

            _fmpz_add_fast(MAG_EXPREF(z), e, 1);

            if (shift >= MAG_BITS)
                MAG_MAN(z) = MAG_ONE_HALF + LIMB_ONE;
            else
                MAG_MAN(z) = MAG_ONE_HALF + (MAG_MAN(x) >> (shift + 1)) + LIMB_ONE;
        }

        MAG_ADJUST_ONE_TOO_LARGE(z);
    }
}

