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
mag_set_d_2exp_fmpz(mag_t z, double c, const fmpz_t exp)
{
    if (c == 0.0)
    {
        mag_zero(z);
    }
    else if (c > 1e300 || c < 0.0) /* not implemented */
    {
        printf("mag_set_d_2exp_fmpz\n");
        abort();
    }
    else
    {
        int cexp, fix;
        mp_limb_t man;

        c = frexp(c, &cexp);

        man = (mp_limb_t)(c * (double)(LIMB_ONE << MAG_BITS)) + 1;

        fix = man >> (MAG_BITS);
        man = (man >> fix) + fix;  /* XXX: need +fix? */
        MAG_MAN(z) = man;
        _fmpz_add_fast(MAG_EXPREF(z), exp, cexp + fix);
    }
}
