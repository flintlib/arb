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
#include "double_extras.h"

void
mag_log1p(mag_t z, const mag_t x)
{
    if (mag_is_special(x))
    {
        if (mag_is_zero(x))
            mag_zero(z);
        else
            mag_inf(z);
    }
    else
    {
        fmpz exp = MAG_EXP(x);

        if (!COEFF_IS_MPZ(exp))
        {
            /* Quick bound by x */
            if (exp < -10)
            {
                mag_set(z, x);
                return;
            }
            else if (exp < 1000)
            {
                double t;
                t = ldexp(MAG_MAN(x), exp - MAG_BITS);
                t = (1.0 + t) * (1 + 1e-14);
                t = mag_d_log_upper_bound(t);
                mag_set_d(z, t);
                return;
            }
        }
        else if (fmpz_sgn(MAG_EXPREF(x)) < 0)
        {
            /* Quick bound by x */
            mag_set(z, x);
            return;
        }

        /* Now we must have x >= 2^1000 */
        /* Use log(2^(exp-1) * (2*v)) = exp*log(2) + log(2*v) */
        {
            double t;
            fmpz_t b;
            mag_t u;

            mag_init(u);
            fmpz_init(b);

            /* incrementing the mantissa gives an upper bound for x+1 */
            t = ldexp(MAG_MAN(x) + 1, 1 - MAG_BITS);
            t = mag_d_log_upper_bound(t);
            mag_set_d(u, t);

            /* log(2) < 744261118/2^30 */
            _fmpz_add_fast(b, MAG_EXPREF(x), -1);
            fmpz_mul_ui(b, b, 744261118);
            mag_set_fmpz(z, b);
            _fmpz_add_fast(MAG_EXPREF(z), MAG_EXPREF(z), -30);

            mag_add(z, z, u);

            mag_clear(u);
            fmpz_clear(b);
        }
    }
}

