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

    Copyright (C) 2015 Fredrik Johansson

******************************************************************************/

#include "mag.h"

void
mag_geom_series(mag_t res, const mag_t x, ulong n)
{
    if (mag_is_zero(x))
    {
        if (n == 0)
            mag_one(res);
        else
            mag_zero(res);
    }
    else if (mag_is_inf(x))
    {
        mag_inf(res);
    }
    else
    {
        mag_t t;
        mag_init(t);
        mag_one(t);
        mag_sub_lower(t, t, x);

        if (mag_is_zero(t))
        {
            mag_inf(res);
        }
        else
        {
            mag_pow_ui(res, x, n);
            mag_div(res, res, t);
        }

        mag_clear(t);
    }
}

