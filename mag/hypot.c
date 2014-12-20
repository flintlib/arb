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
mag_hypot(mag_t z, const mag_t x, const mag_t y)
{
    if (mag_is_zero(y))
    {
        mag_set(z, x);
    }
    else if (mag_is_zero(x))
    {
        mag_set(z, y);
    }
    else
    {
        mag_t t;
        mag_init(t);
        mag_mul(t, x, x);
        mag_addmul(t, y, y);
        mag_sqrt(z, t);
        mag_clear(t);
    }
}

