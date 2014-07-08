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
mag_sqrt(mag_t y, const mag_t x)
{
    if (mag_is_special(x))
    {
        mag_set(y, x);
    }
    else
    {
        double t;
        fmpz e;

        t = MAG_MAN(x) * ldexp(1.0, -MAG_BITS);
        e = MAG_EXP(x);

        if (!COEFF_IS_MPZ(e))
        {
            if (e % 2 != 0)
            {
                e = (e - 1) >> 1;
                t *= 2.0;
            }
            else
            {
                e >>= 1;
            }
            t = sqrt(t) * (1 + 1e-13);
            mag_set_d_2exp_fmpz(y, t, &e);
        }
        else
        {
            if (fmpz_is_odd(MAG_EXPREF(x)))
                t *= 2.0;
            fmpz_fdiv_q_2exp(MAG_EXPREF(y), MAG_EXPREF(x), 1);
            t = sqrt(t) * (1 + 1e-13);
            mag_set_d_2exp_fmpz(y, t, MAG_EXPREF(y));
        }
    }
}
