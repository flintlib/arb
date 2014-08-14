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

static const double inverse_factorials[] = {
    1.0,
    1.0,
    0.5,
    0.16666666666666666667,
    0.041666666666666666667,
    0.0083333333333333333333,
    0.0013888888888888888889,
    0.0001984126984126984127,
    0.000024801587301587301587,
    2.7557319223985890653e-6,
    2.7557319223985890653e-7,
    2.5052108385441718775e-8,
    2.0876756987868098979e-9,
    1.6059043836821614599e-10,
    1.1470745597729724714e-11,
    7.6471637318198164759e-13
};

static __inline__ double
_mag_d_exp_upper_reduced(double u)
{
    if (u < -0.375 || u > 0.375)
        abort();

    return d_polyval(inverse_factorials, 11, u) + 1e-12;
}

void
mag_exp_maglim(mag_t y, const mag_t x, long maglim)
{
    if (mag_is_special(x))
    {
        if (mag_is_zero(x))
            mag_one(y);
        else
            mag_inf(y);
    }
    else if (COEFF_IS_MPZ(MAG_EXP(x)))
    {
        if (fmpz_sgn(MAG_EXPREF(x)) > 0)
        {
            mag_inf(y);
        }
        else
        {
            MAG_MAN(y) = MAG_ONE_HALF + 1;
            fmpz_one(MAG_EXPREF(y));
        }
    }
    else
    {
        long e = MAG_EXP(x);

        if (e <= -MAG_BITS) /* assumes MAG_BITS == 30 */
        {
            MAG_MAN(y) = MAG_ONE_HALF + 1;
            fmpz_one(MAG_EXPREF(y));
        }
        else if (e <= -(MAG_BITS / 2))  /* assumes MAG_BITS == 30 */
        {
            MAG_MAN(y) = MAG_ONE_HALF + (MAG_MAN(x) >> (1-e)) + 2;
            fmpz_one(MAG_EXPREF(y));
        }
        else if (e < 24)
        {
            double t, u;
            ulong n;

            t = ldexp(MAG_MAN(x), e - MAG_BITS);

            /* does not need to be exact */
            n = (ulong)(t * 1.4426950408889634074 + 0.5);
            /* here u must be rounded up */
            u = t - n * (0.69314718055994530942 * (1.0 - 1e-13)) + 1e-13;

            u = _mag_d_exp_upper_reduced(u);
            fmpz_set_ui(MAG_EXPREF(y), n);
            mag_set_d_2exp_fmpz(y, u, MAG_EXPREF(y));
        }
        else if (e > maglim)
        {
            mag_inf(y);
        }
        else
        {
            /* we really want a multiprecision algorithm
               here for huge n, but for most purposes, it's fine to just
               get a few leading digits of the *exponent* accurately */
            fmpz_t t;
            fmpz_init(t);

            fmpz_set_ui(t, MAG_MAN(x));

            if (e >= MAG_BITS)
                fmpz_mul_2exp(t, t, e - MAG_BITS);
            else
                fmpz_cdiv_q_2exp(t, t, MAG_BITS - e);

            /* upper bound for e */
            MAG_MAN(y) = 729683223;
            fmpz_set_ui(MAG_EXPREF(y), 2);

            mag_pow_fmpz(y, y, t);
            fmpz_clear(t);
        }
    }
}

