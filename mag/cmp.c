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

int
mag_cmp(const mag_t x, const mag_t y)
{
    int c;

    if (mag_equal(x, y))
        return 0;

    if (mag_is_special(x) || mag_is_special(y))
    {
        if (mag_is_zero(x)) return -1;
        if (mag_is_zero(y)) return 1;
        if (mag_is_inf(x)) return 1;
        if (mag_is_inf(y)) return -1;
    }

    c = fmpz_cmp(MAG_EXPREF(x), MAG_EXPREF(y));

    if (c == 0)
        return (MAG_MAN(x) < MAG_MAN(y)) ? -1 : 1;

    return (c < 0) ? -1 : 1;
}

