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

#define BIG_EXPONENT_BITS 20
#define BIG_EXPONENT (1L << BIG_EXPONENT_BITS)

void
fmpr_log1p_ubound(fmpr_t y, const fmpr_t x)
{
    if (fmpr_is_special(x))
    {
        if (fmpr_is_zero(x))
            fmpr_zero(y);
        else if (fmpr_is_pos_inf(x))
            fmpr_pos_inf(y);
        else
            fmpr_nan(y);
    }
    else if (fmpr_sgn(x) < 0)
    {
        fmpr_nan(y);
    }
    else
    {
        fmpz_t exp;
        fmpz_init(exp);
        fmpz_add_ui(exp, fmpr_expref(x), fmpz_bits(fmpr_manref(x)));

        if (fmpz_cmp_si(exp, -8) < 0)
        {
            fmpr_set(y, x);
        }
        else if (fmpz_cmp_si(exp, BIG_EXPONENT) > 0)
        {
            /* bound 1+x by 2^n */
            fmpr_add_ui(y, x, 1, FMPRB_RAD_PREC, FMPR_RND_UP);
            fmpz_add_ui(fmpr_manref(y), fmpr_expref(y), fmpz_bits(fmpr_manref(y)));

            /* log(2) < 11629080 / 2^24 */
            fmpz_mul_ui(fmpr_manref(y), fmpr_manref(y), 11629080);
            fmpz_set_si(fmpr_expref(y), -24);
            _fmpr_normalise(fmpr_manref(y), fmpr_expref(y), FMPRB_RAD_PREC, FMPR_RND_UP);
        }
        else
        {
            fmpr_log1p(y, x, FMPRB_RAD_PREC, FMPR_RND_UP);
        }

        fmpz_clear(exp);
    }
}

void
fmprb_log_fmpr(fmprb_t y, const fmpr_t x, long prec)
{
    long r;
    fmpz_t exp;

    if (fmpr_is_special(x))
    {
        if (fmpr_is_pos_inf(x))
        {
            fmpr_pos_inf(fmprb_midref(y));
            fmpr_zero(fmprb_radref(y));
        }
        else
        {
            fmpr_nan(fmprb_midref(y));
            fmpr_pos_inf(fmprb_radref(y));
        }
        return;
    }

    if (fmpz_sgn(fmpr_manref(x)) < 0)
    {
        fmpr_nan(fmprb_midref(y));
        fmpr_pos_inf(fmprb_radref(y));
        return;
    }

    if (fmpz_is_one(fmpr_manref(x)))
    {
        fmpz_init(exp);
        fmpz_set(exp, fmpr_expref(x));
        fmprb_const_log2(y, prec + 2);
        fmprb_mul_fmpz(y, y, exp, prec);
        fmpz_clear(exp);
        return;
    }

    fmpz_init(exp);
    fmpz_add_ui(exp, fmpr_expref(x), fmpz_bits(fmpr_manref(x)));

    if (*exp >= -BIG_EXPONENT && *exp <= BIG_EXPONENT)
    {
        r = fmpr_log(fmprb_midref(y), x, prec, FMPR_RND_DOWN);
        fmpr_set_error_result(fmprb_radref(y), fmprb_midref(y), r);
    }
    else
    {
        fmpr_t t;
        fmprb_t c;
        long wp;

        fmpr_init(t);
        fmprb_init(c);

        fmpz_neg(exp, exp);
        fmpr_mul_2exp_fmpz(t, x, exp);

        wp = prec + 4 - fmpz_bits(exp);
        wp = FLINT_MAX(wp, 4);

        r = fmpr_log(fmprb_midref(y), t, wp, FMPR_RND_DOWN);
        fmpr_set_error_result(fmprb_radref(y), fmprb_midref(y), r);

        fmprb_const_log2(c, prec + 4);
        fmprb_mul_fmpz(c, c, exp, prec + 4);

        fmprb_sub(y, y, c, prec);

        fmpr_clear(t);
        fmprb_clear(c);
    }

    fmpz_clear(exp);
}

void
fmprb_log(fmprb_t y, const fmprb_t x, long prec)
{
    if (fmprb_is_exact(x))
    {
        fmprb_log_fmpr(y, fmprb_midref(x), prec);
    }
    else
    {
        /*
        Let the input be [a-b, a+b]. We require a > b >= 0 (otherwise the
        interval contains zero or a negative number and the logarithm is not
        defined). The error is largest at a-b, and we have

        log(a) - log(a-b) = log(1 + b/(a-b)).
        */
        fmpr_t err;
        fmpr_init(err);
        fmpr_sub(err, fmprb_midref(x), fmprb_radref(x), FMPRB_RAD_PREC, FMPR_RND_DOWN);

        if (fmpr_sgn(err) <= 0)
        {
            fmpr_pos_inf(err);
        }
        else
        {
            fmpr_div(err, fmprb_radref(x), err, FMPRB_RAD_PREC, FMPR_RND_UP);
            fmpr_log1p_ubound(err, err);
        }

        fmprb_log_fmpr(y, fmprb_midref(x), prec);
        fmpr_add(fmprb_radref(y), fmprb_radref(y), err, FMPRB_RAD_PREC, FMPR_RND_UP);
        fmpr_clear(err);
    }

    fmprb_adjust(y);
}

void
fmprb_log_ui(fmprb_t z, ulong x, long prec)
{
    fmpr_t t;
    fmpr_init(t);
    fmpr_set_ui(t, x);
    fmprb_log_fmpr(z, t, prec);
    fmpr_clear(t);
}

void
fmprb_log_fmpz(fmprb_t z, const fmpz_t x, long prec)
{
    fmpr_t t;
    fmpr_init(t);
    fmpr_set_fmpz(t, x);
    fmprb_log_fmpr(z, t, prec);
    fmpr_clear(t);
}
