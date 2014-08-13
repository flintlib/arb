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
mag_randtest_special(mag_t x, flint_rand_t state, long expbits)
{
    switch (n_randint(state, 32))
    {
        case 0:
            mag_zero(x);
            break;
        case 1:
            mag_inf(x);
            break;
        case 2:
            MAG_MAN(x) = (LIMB_ONE << MAG_BITS) - 1;
            fmpz_randtest(MAG_EXPREF(x), state, expbits);
            break;
        case 3:
            MAG_MAN(x) = LIMB_ONE << (MAG_BITS - 1);
            fmpz_randtest(MAG_EXPREF(x), state, expbits);
            break;
        default:
            MAG_MAN(x) = (n_randtest(state) >> (FLINT_BITS - MAG_BITS));
            MAG_MAN(x) |= (LIMB_ONE << (MAG_BITS - 1));
            fmpz_randtest(MAG_EXPREF(x), state, expbits);
    }
}

void
mag_randtest(mag_t x, flint_rand_t state, long expbits)
{
    mag_randtest_special(x, state, expbits);
    if (mag_is_inf(x))
        mag_zero(x);
}

