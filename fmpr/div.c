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

#include "fmpr.h"

static void
_fmpr_div_special(fmpr_t z, const fmpr_t x, const fmpr_t y)
{
    if ((fmpr_is_zero(x) && !fmpr_is_zero(y) && !fmpr_is_nan(y)) ||
        (fmpr_is_inf(y) && !fmpr_is_special(x)))
    {
        fmpr_zero(z);
    }
    else if (fmpr_is_zero(y) || (fmpr_is_special(x) && fmpr_is_special(y)) ||
        fmpr_is_nan(x) || fmpr_is_nan(y))
    {
        fmpr_nan(z);
    }
    else if (fmpr_sgn(x) == fmpr_sgn(y))
        fmpr_pos_inf(z);
    else
        fmpr_neg_inf(z);
}

void
fmpr_div(fmpr_t z, const fmpr_t x, const fmpr_t y, long prec, fmpr_rnd_t rnd)
{
    if (fmpr_is_special(x) || fmpr_is_special(y))
    {
        _fmpr_div_special(z, x, y);
        return;
    }

    /* division by power of two <=> shift exponents */
    if (fmpz_is_pm1(fmpr_manref(y)))
    {
        if (fmpz_is_one(fmpr_manref(y)))
            fmpz_set(fmpr_manref(z), fmpr_manref(x));
        else
            fmpz_neg(fmpr_manref(z), fmpr_manref(x));
        fmpz_sub(fmpr_expref(z), fmpr_expref(x), fmpr_expref(y));
        _fmpr_normalise(fmpr_manref(z), fmpr_expref(z), prec, rnd);
        return;
    }
    else
    {
        long extra;
        int negative;
        fmpz_t t, rem;

        /* todo: work out exact needed shift */
        extra = prec - fmpz_bits(fmpr_manref(x)) + fmpz_bits(fmpr_manref(y)) + 5;
        extra = FLINT_MAX(extra, 5);

        negative = fmpz_sgn(fmpr_manref(x)) != fmpz_sgn(fmpr_manref(y));

        fmpz_init(rem);
        fmpz_init(t);

        fmpz_mul_2exp(t, fmpr_manref(x), extra);
        fmpz_tdiv_qr(fmpr_manref(z), rem, t, fmpr_manref(y));

        if (!fmpz_is_zero(rem))
        {
            fmpz_mul_2exp(fmpr_manref(z), fmpr_manref(z), 1);

            if (negative)
                fmpz_sub_ui(fmpr_manref(z), fmpr_manref(z), 1);
            else
                fmpz_add_ui(fmpr_manref(z), fmpr_manref(z), 1);

            extra++;
        }

        fmpz_clear(rem);
        fmpz_clear(t);

        fmpz_sub(fmpr_expref(z), fmpr_expref(x), fmpr_expref(y));
        fmpz_sub_ui(fmpr_expref(z), fmpr_expref(z), extra);

        _fmpr_normalise(fmpr_manref(z), fmpr_expref(z), prec, rnd);
    }
}

