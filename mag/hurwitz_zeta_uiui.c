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

    Copyright (C) 2016 Fredrik Johansson

******************************************************************************/

#include "mag.h"

void
mag_hurwitz_zeta_uiui(mag_t res, ulong s, ulong a)
{
    if (s <= 1 || a == 0)
    {
        mag_inf(res);
    }
    else
    {
        mag_t t;
        mag_init(t);
        mag_one(t);

        /* 1/a + 1/(a^(s-1)*(s-1)) */
        mag_set_ui_lower(res, a);
        mag_pow_ui_lower(res, res, s - 1);
        mag_mul_ui_lower(res, res, s - 1);
        mag_div(res, t, res);

        mag_div_ui(t, t, a);
        mag_add(res, res, t);

        mag_clear(t);
    }
}

