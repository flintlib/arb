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
mag_root(mag_t y, const mag_t x, ulong n)
{
    if (n == 0)
    {
        mag_inf(y);
    }
    else if (n == 1 || mag_is_special(x))
    {
        mag_set(y, x);
    }
    else if (n == 2)
    {
        mag_sqrt(y, x);
    }
    else if (n == 4)
    {
        mag_sqrt(y, x);
        mag_sqrt(y, y);
    }
    else
    {
        fmpz_t e, f;

        fmpz_init_set_ui(e, MAG_BITS);
        fmpz_init(f);

        fmpz_sub(e, e, MAG_EXPREF(x));
        fmpz_cdiv_q_ui(e, e, n);
        fmpz_mul_ui(f, e, n);
        mag_mul_2exp_fmpz(y, x, f);
        mag_log1p(y, y);
        mag_div_ui(y, y, n);
        mag_exp(y, y);
        fmpz_neg(e, e);
        mag_mul_2exp_fmpz(y, y, e);

        fmpz_clear(e);
        fmpz_clear(f);
    }
}

