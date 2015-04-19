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

#include "arf.h"

void
arf_get_mag(mag_t y, const arf_t x)
{
    if (arf_is_zero(x))
    {
        mag_zero(y);
    }
    else if (arf_is_special(x))
    {
        mag_inf(y);
    }
    else
    {
        mp_limb_t t, u;

        ARF_GET_TOP_LIMB(t, x);
        t = (t >> (FLINT_BITS - MAG_BITS)) + LIMB_ONE;
        u = t >> MAG_BITS;
        t = (t >> u) + u;
        _fmpz_add_fast(MAG_EXPREF(y), ARF_EXPREF(x), u);
        MAG_MAN(y) = t;
    }
}

void
arf_get_mag_lower(mag_t y, const arf_t x)
{
    if (arf_is_zero(x))
    {
        mag_zero(y);
    }
    else if (arf_is_special(x))
    {
        if (arf_is_nan(x))
            mag_zero(y);
        else
            mag_inf(y);
    }
    else
    {
        mp_limb_t t;
        ARF_GET_TOP_LIMB(t, x);
        MAG_MAN(y) = t >> (FLINT_BITS - MAG_BITS);
        _fmpz_set_fast(MAG_EXPREF(y), ARF_EXPREF(x));
    }
}

