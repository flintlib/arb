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

static long _fmpr_add_special(fmpr_t z, const fmpr_t x, const fmpr_t y, long prec, fmpr_rnd_t rnd)
{
    if (fmpr_is_zero(x))
    {
        if (fmpr_is_zero(y))
        {
            fmpr_zero(z);
            return FMPR_RESULT_EXACT;
        }
        else
            return fmpr_set_round(z, y, prec, rnd);
    }
    else if (fmpr_is_zero(y))
    {
        return fmpr_set_round(z, x, prec, rnd);
    }
    else if (fmpr_is_nan(x) || fmpr_is_nan(y)
        || (fmpr_is_pos_inf(x) && fmpr_is_neg_inf(y))
        || (fmpr_is_neg_inf(x) && fmpr_is_pos_inf(y)))
    {
        fmpr_nan(z);
        return FMPR_RESULT_EXACT;
    }
    else if (fmpr_is_special(x))
    {
        fmpr_set(z, x);
        return FMPR_RESULT_EXACT;
    }
    else
    {
        fmpr_set(z, y);
        return FMPR_RESULT_EXACT;
    }
}

long
fmpr_add_naive(fmpr_t z, const fmpr_t x, const fmpr_t y, long prec, fmpr_rnd_t rnd)
{
    long shift, xsize, ysize;

    if (fmpr_is_special(x) || fmpr_is_special(y))
    {
        return _fmpr_add_special(z, x, y, prec, rnd);
    }

    shift = _fmpz_sub_small(fmpr_expref(x), fmpr_expref(y));

    if (shift == 0)
    {
        fmpz_add(fmpr_manref(z), fmpr_manref(x), fmpr_manref(y));
        fmpz_set(fmpr_expref(z), fmpr_expref(x));
    }
    else if (shift > 0)
    {
        ysize = _fmpz_size(fmpr_manref(y)) * FLINT_BITS;

        /* x and y do not overlap */
        if (shift > ysize && prec != FMPR_PREC_EXACT)
        {
            /* y does not overlap with result */
            if (ysize + prec - (long) fmpz_bits(fmpr_manref(x)) < shift)
            {
                return _fmpr_add_eps(z, x, fmpz_sgn(fmpr_manref(y)), prec, rnd);
            }
        }

        fmpz_add_mul2exp(fmpr_manref(z), fmpr_manref(y), fmpr_manref(x), shift);
        fmpz_set(fmpr_expref(z), fmpr_expref(y));
    }
    else
    {
        shift = -shift;

        xsize = _fmpz_size(fmpr_manref(x)) * FLINT_BITS;

        /* x and y do not overlap */
        if (shift > xsize && prec != FMPR_PREC_EXACT)
        {
            /* y does not overlap with result */
            if (xsize + prec - (long) fmpz_bits(fmpr_manref(y)) < shift)
            {
                return _fmpr_add_eps(z, y, fmpz_sgn(fmpr_manref(x)), prec, rnd);
            }
        }

        fmpz_add_mul2exp(fmpr_manref(z), fmpr_manref(x), fmpr_manref(y), shift);
        fmpz_set(fmpr_expref(z), fmpr_expref(x));
    }

    return _fmpr_normalise(fmpr_manref(z), fmpr_expref(z), prec, rnd);
}

