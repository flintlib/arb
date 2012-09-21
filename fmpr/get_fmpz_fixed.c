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

int
fmpr_get_fmpz_fixed_fmpz(fmpz_t y, const fmpr_t x, const fmpz_t exp)
{
    long shift;

    if (fmpr_is_zero(x))
    {
        fmpz_zero(y);
        return 0;
    }

    shift = _fmpz_sub_small(fmpr_expref(x), exp);

    if (shift >= 0)
    {
        fmpz_mul_2exp(y, fmpr_manref(x), shift);
        return 0;
    }
    else
    {
        fmpz_tdiv_q_2exp(y, fmpr_manref(x), -shift);
        return 1;
    }
}

int
fmpr_get_fmpz_fixed_si(fmpz_t y, const fmpr_t x, long exp)
{
    int result;
    fmpz_t t;
    fmpz_init(t);
    fmpz_set_si(t, exp);
    result = fmpr_get_fmpz_fixed_fmpz(y, x, t);
    fmpz_clear(t);
    return result;
}
