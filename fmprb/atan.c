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

static long
_fmpr_atan(fmpr_t y, const fmpr_t x, long prec, fmpr_rnd_t rnd)
{
    long r;
    CALL_MPFR_FUNC(r, mpfr_atan, y, x, prec, rnd);
    return r;
}

void
fmprb_atan_fmpr(fmprb_t z, const fmpr_t x, long prec)
{
    if (fmpr_is_special(x))
    {
        if (fmpr_is_zero(x))
        {
            fmprb_zero(z);
        }
        else if (fmpr_is_nan(x))
        {
            fmpr_nan(fmprb_midref(z));
            fmpr_pos_inf(fmprb_radref(z));
        }
        else if (fmpr_is_pos_inf(x))
        {
            fmprb_const_pi(z, prec);
            fmprb_mul_2exp_si(z, z, -1);
        }
        else if (fmpr_is_neg_inf(x))
        {
            fmprb_const_pi(z, prec);
            fmprb_mul_2exp_si(z, z, -1);
            fmprb_neg(z, z);
        }
    }
    else
    {
        long r;
        fmpz_t mag;

        fmpz_init(mag);

        /* 2^(mag-1) <= |x| < 2^mag */
        fmpz_add_ui(mag, fmpr_expref(x), fmpz_bits(fmpr_manref(x)));

        /* atan(x) = x + eps, |eps| < x^3 */
        if (fmpz_cmp_si(mag, -(prec/3) - 2) < 0)
        {
            fmpz_mul_ui(mag, mag, 3);
            fmprb_set_fmpr(z, x);
            fmprb_set_round(z, z, prec);
            fmprb_add_error_2exp_fmpz(z, mag);
        }
        /* atan(x) = pi/2 - eps, eps < 1/x <= 2^(1-mag) */
        /* TODO: also use atan(x) = pi/2 - 1/x + eps, eps < 1/x^3 */
        else if (fmpz_cmp_si(mag, prec + 2) > 0)
        {
            fmpz_neg(mag, mag);
            fmpz_add_ui(mag, mag, 1);
            if (fmpr_sgn(x) > 0)
            {
                fmprb_const_pi(z, prec);
            }
            else
            {
                fmprb_const_pi(z, prec);
                fmprb_neg(z, z);
            }
            fmprb_mul_2exp_si(z, z, -1);
            fmprb_add_error_2exp_fmpz(z, mag);
        }
        else
        {
            r = _fmpr_atan(fmprb_midref(z), x, prec, FMPR_RND_DOWN);
            fmpr_set_error_result(fmprb_radref(z), fmprb_midref(z), r);
        }

        fmpz_clear(mag);
    }
}

void
fmprb_atan(fmprb_t z, const fmprb_t x, long prec)
{
    if (fmprb_is_exact(x))
    {
        fmprb_atan_fmpr(z, fmprb_midref(x), prec);
    }
    else
    {
        fmpr_t t, u;

        fmpr_init(t);
        fmpr_init(u);

        if (fmpr_sgn(fmprb_midref(x)) >= 0)
        {
            fmpr_sub(t, fmprb_midref(x), fmprb_radref(x), FMPRB_RAD_PREC, FMPR_RND_DOWN);
        }
        else
        {
            fmpr_add(t, fmprb_midref(x), fmprb_radref(x), FMPRB_RAD_PREC, FMPR_RND_DOWN);
            fmpr_neg(t, t);
        }

        if (fmpr_sgn(t) > 0)
        {
            fmpr_mul(t, t, t, FMPRB_RAD_PREC, FMPR_RND_DOWN);
            fmpr_add_ui(t, t, 1UL, FMPRB_RAD_PREC, FMPR_RND_DOWN);
            fmpr_div(t, fmprb_radref(x), t, FMPRB_RAD_PREC, FMPR_RND_UP);
        }
        else
        {
            fmpr_set(t, fmprb_radref(x));
        }

        fmprb_atan_fmpr(z, fmprb_midref(x), prec);
        fmpr_add(fmprb_radref(z), fmprb_radref(z), t, FMPRB_RAD_PREC, FMPR_RND_UP);

        fmpr_clear(t);
        fmpr_clear(u);
    }

    fmprb_adjust(z);
}

