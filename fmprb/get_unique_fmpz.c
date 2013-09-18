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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "fmprb.h"

int
fmprb_get_unique_fmpz(fmpz_t z, const fmprb_t x)
{
    if (!fmprb_is_finite(x))
    {
        return 0;
    }
    else if (fmpr_is_zero(fmprb_radref(x)))
    {
        /* x = b*2^e, e >= 0 */
        if (fmpr_is_int(fmprb_midref(x)))
        {
            if (!fmpz_fits_si(fmpr_expref(fmprb_midref(x))))
            {
                printf("fmprb_get_unique_fmpz: too large shift\n");
                abort();
            }

            fmpz_mul_2exp(z, fmpr_manref(fmprb_midref(x)),
                fmpz_get_si(fmpr_expref(fmprb_midref(x))));

            return 1;
        }
        else
        {
            return 0;
        }
    }
    /* if the radius is >= 1, there are at least two integers */
    else if (fmpr_cmp_2exp_si(fmprb_radref(x), 0) >= 0)
    {
        return 0;
    }
    /* there are 0 or 1 integers if the radius is < 1 */
    else
    {
        /* if the midpoint is exactly an integer, it is what we want */
        if (fmpr_is_int(fmprb_midref(x)))
        {
            if (!fmpz_fits_si(fmpr_expref(fmprb_midref(x))))
            {
                printf("fmprb_get_unique_fmpz: too large shift\n");
                abort();
            }

            fmpz_mul_2exp(z, fmpr_manref(fmprb_midref(x)),
                fmpz_get_si(fmpr_expref(fmprb_midref(x))));

            return 1;
        }
        /* if the radius is tiny, it can't be an integer */
        else if (!fmpz_fits_si(fmpr_expref(fmprb_radref(x))))
        {
            return 0;
        }
        else
        {
            fmpz_t a, b, exp;
            int res;

            fmpz_init(a);
            fmpz_init(b);
            fmpz_init(exp);

            fmprb_get_interval_fmpz_2exp(a, b, exp, x);

            fmpz_cdiv_q_2exp(a, a, -fmpz_get_si(exp));
            fmpz_fdiv_q_2exp(b, b, -fmpz_get_si(exp));

            res = fmpz_equal(a, b);

            if (res)
                fmpz_set(z, a);

            fmpz_clear(a);
            fmpz_clear(b);
            fmpz_clear(exp);

            return res;
        }
    }
}

