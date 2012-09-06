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
_fmpr_mul_special(fmpr_t z, const fmpr_t x, const fmpr_t y)
{
    if (fmpr_is_zero(x))
    {
        if (!fmpr_is_special(y) || fmpr_is_zero(y))
            fmpr_zero(z);
        else
            fmpr_nan(z);
        return;
    }

    if (fmpr_is_zero(y))
    {
        if (!fmpr_is_special(x))
            fmpr_zero(z);
        else
            fmpr_nan(z);
        return;
    }

    if ((fmpr_is_inf(x) && (fmpr_is_inf(y) || !fmpr_is_special(y))) ||
        (fmpr_is_inf(y) && !fmpr_is_special(x)))
    {
        if (fmpr_sgn(x) == fmpr_sgn(y))
            fmpr_pos_inf(z);
        else
            fmpr_neg_inf(z);
        return;
    }

    fmpr_nan(z);
}

long
fmpr_mul(fmpr_t z, const fmpr_t x, const fmpr_t y, long prec, fmpr_rnd_t rnd)
{
    if (fmpr_is_special(x) || fmpr_is_special(y))
    {
        _fmpr_mul_special(z, x, y);
        return FMPR_RESULT_EXACT;
    }

    fmpz_mul(fmpr_manref(z), fmpr_manref(x), fmpr_manref(y));
    fmpz_add(fmpr_expref(z), fmpr_expref(x), fmpr_expref(y));
    return _fmpr_normalise(fmpr_manref(z), fmpr_expref(z), prec, rnd);
}

long
fmpr_mul_ui(fmpr_t z, const fmpr_t x, ulong y, long prec, fmpr_rnd_t rnd)
{
    fmpr_t t; long r;
    fmpr_init(t);
    fmpr_set_ui(t, y);
    r = fmpr_mul(z, x, t, prec, rnd);
    fmpr_clear(t);
    return r;
}

long
fmpr_mul_si(fmpr_t z, const fmpr_t x, long y, long prec, fmpr_rnd_t rnd)
{
    fmpr_t t; long r;
    fmpr_init(t);
    fmpr_set_si(t, y);
    r = fmpr_mul(z, x, t, prec, rnd);
    fmpr_clear(t);
    return r;
}

long
fmpr_mul_fmpz(fmpr_t z, const fmpr_t x, const fmpz_t y, long prec, fmpr_rnd_t rnd)
{
    fmpr_t t; long r;
    fmpr_init(t);
    fmpr_set_fmpz(t, y);
    r = fmpr_mul(z, x, t, prec, rnd);
    fmpr_clear(t);
    return r;
}
