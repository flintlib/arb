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

long
fmpr_sqrt(fmpr_t y, const fmpr_t x, long prec, fmpr_rnd_t rnd)
{
    long shift, man_shift, r, bc;
    fmpz_t rem;

    if (fmpr_is_special(x))
    {
        if (fmpr_is_zero(x))
            fmpr_zero(y);
        else if (fmpr_is_pos_inf(x))
            fmpr_pos_inf(y);
        else
            fmpr_nan(y);

        return FMPR_RESULT_EXACT;
    }

    if (fmpr_sgn(x) < 0)
    {
        fmpr_nan(y);
        return FMPR_RESULT_EXACT;
    }

    /* special case: 4^n */
    /* TODO: process all small exact square roots efficiently */
    if (fmpz_is_one(fmpr_manref(x)) && fmpz_is_even(fmpr_expref(x)))
    {
        r = fmpr_set_round(y, x, prec, rnd);
        fmpz_tdiv_q_2exp(fmpr_expref(y), fmpr_expref(y), 1);
        return r;
    }

    fmpz_set(fmpr_expref(y), fmpr_expref(x));

    man_shift = 0;

    bc = fmpz_bits(fmpr_manref(x));

    if (fmpz_is_odd(fmpr_expref(y)))
    {
        fmpz_sub_ui(fmpr_expref(y), fmpr_expref(y), 1UL);
        man_shift += 1;
        bc += 1;
    }

    shift = FLINT_MAX(4, 2 * prec - bc + 4);
    shift += shift & 1;

    fmpz_init(rem);
    fmpz_mul_2exp(fmpr_manref(y), fmpr_manref(x), shift + man_shift);
    fmpz_sqrtrem(fmpr_manref(y), rem, fmpr_manref(y));

    if (!fmpz_is_zero(rem) && rnd != FMPR_RND_FLOOR && rnd != FMPR_RND_DOWN)
    {
        fmpz_mul_2exp(fmpr_manref(y), fmpr_manref(y), 1);
        fmpz_add_ui(fmpr_manref(y), fmpr_manref(y), 1UL);
        shift += 2;
    }

    fmpz_clear(rem);

    fmpz_sub_ui(fmpr_expref(y), fmpr_expref(y), shift);
    fmpz_tdiv_q_2exp(fmpr_expref(y), fmpr_expref(y), 1);

    return _fmpr_normalise(fmpr_manref(y), fmpr_expref(y), prec, rnd);
}


long
fmpr_sqrt_ui(fmpr_t z, ulong x, long prec, fmpr_rnd_t rnd)
{
    fmpr_t t; long r;
    fmpr_init(t);
    fmpr_set_ui(t, x);
    r = fmpr_sqrt(z, t, prec, rnd);
    fmpr_clear(t);
    return r;
}

long
fmpr_sqrt_fmpz(fmpr_t z, const fmpz_t x, long prec, fmpr_rnd_t rnd)
{
    fmpr_t t; long r;
    fmpr_init(t);
    fmpr_set_fmpz(t, x);
    r = fmpr_sqrt(z, t, prec, rnd);
    fmpr_clear(t);
    return r;
}
