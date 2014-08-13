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
mag_pow_fmpz(mag_t z, const mag_t x, const fmpz_t e)
{
    if (fmpz_sgn(e) < 0)
    {
        abort();
    }
    else if (!COEFF_IS_MPZ(*e))
    {
        mag_pow_ui(z, x, fmpz_get_ui(e));
    }
    else
    {
        mag_t y;
        mp_srcptr elimbs;
        long i, bits;

        mag_init_set(y, x);
        bits = fmpz_bits(e);
        elimbs = COEFF_TO_PTR(*e)->_mp_d;

        for (i = bits - 2; i >= 0; i--)
        {
            mag_mul(y, y, y);

            if ((elimbs[i / FLINT_BITS] >> (i % FLINT_BITS)) & 1)
                mag_mul(y, y, x);
        }

        mag_swap(z, y);
        mag_clear(y);
    }
}

