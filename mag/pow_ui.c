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
mag_pow_ui(mag_t z, const mag_t x, ulong e)
{
    if (mag_is_inf(x))
    {
        mag_inf(z);
    }
    else if (e <= 2)
    {
        if (e == 0)
            mag_one(z);
        else if (e == 1)
            mag_set(z, x);
        else
            mag_mul(z, x, x);
    }
    else
    {
        mag_t y;
        int i, bits;

        mag_init_set(y, x);

        bits = FLINT_BIT_COUNT(e);

        for (i = bits - 2; i >= 0; i--)
        {
            mag_mul(y, y, y);
            if (e & (1UL << i))
                mag_mul(y, y, x);
        }

        mag_swap(z, y);
        mag_clear(y);
    }
}

void
mag_pow_ui_lower(mag_t z, const mag_t x, ulong e)
{
    if (e <= 2)
    {
        if (e == 0)
            mag_one(z);
        else if (e == 1)
            mag_set(z, x);
        else
            mag_mul_lower(z, x, x);
    }
    else if (mag_is_inf(x))
    {
        mag_inf(z);
    }
    else
    {
        mag_t y;
        int i, bits;

        mag_init_set(y, x);

        bits = FLINT_BIT_COUNT(e);

        for (i = bits - 2; i >= 0; i--)
        {
            mag_mul_lower(y, y, y);
            if (e & (1UL << i))
                mag_mul_lower(y, y, x);
        }

        mag_swap(z, y);
        mag_clear(y);
    }
}

