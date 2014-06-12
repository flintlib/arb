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
mag_set_fmpz_2exp_fmpz(mag_t z, const fmpz_t man, const fmpz_t exp)
{
    if (fmpz_is_zero(man))
    {
        mag_zero(z);
    }
    else
    {
        mp_limb_t m;
        long cexp;

        m = fmpz_abs_ubound_ui_2exp(&cexp, man, MAG_BITS);
        MAG_MAN(z) = m;
        _fmpz_add_fast(MAG_EXPREF(z), exp, cexp + MAG_BITS);
    }
}

void
mag_set_fmpz_2exp_fmpz_lower(mag_t z, const fmpz_t man, const fmpz_t exp)
{
    if (fmpz_is_zero(man))
    {
        mag_zero(z);
    }
    else
    {
        mp_limb_t m;
        long cexp;

        m = fmpz_abs_lbound_ui_2exp(&cexp, man, MAG_BITS);
        MAG_MAN(z) = m;
        _fmpz_add_fast(MAG_EXPREF(z), exp, cexp + MAG_BITS);
    }
}

