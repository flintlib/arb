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
mag_expinv(mag_t res, const mag_t x)
{
    if (mag_is_zero(x))
    {
        mag_one(res);
    }
    else if (mag_is_inf(x))
    {
        mag_zero(res);
    }
    else if (fmpz_sgn(MAG_EXPREF(x)) <= 0)
    {
        mag_one(res);
    }
    else if (fmpz_cmp_ui(MAG_EXPREF(x), 2 * MAG_BITS) > 0)
    {
        fmpz_t t;
        fmpz_init(t);

        /* If x > 2^60, exp(-x) < 2^(-2^60 / log(2))  */
        /* -1/log(2) < -369/256 */
        fmpz_set_si(t, -369);
        fmpz_mul_2exp(t, t, 2 * MAG_BITS - 8);

        mag_one(res);
        mag_mul_2exp_fmpz(res, res, t);

        fmpz_clear(t);
    }
    else
    {
        fmpz_t t;
        slong e = MAG_EXP(x);

        fmpz_init(t);
        fmpz_set_ui(t, MAG_MAN(x));

        if (e >= MAG_BITS)
            fmpz_mul_2exp(t, t, e - MAG_BITS);
        else
            fmpz_tdiv_q_2exp(t, t, MAG_BITS - e);

        /* upper bound for 1/e */
        mag_set_ui_2exp_si(res, 395007543, -30);

        mag_pow_fmpz(res, res, t);
        fmpz_clear(t);
    }
}

