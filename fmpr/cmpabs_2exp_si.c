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

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include "fmpr.h"

static __inline__ int sgn(int x)
{
    return (x > 0) - (x < 0);
}

int
fmpr_cmpabs_2exp_si(const fmpr_t x, long e)
{
    long bc;
    int ret;
    fmpz_t t;

    if (fmpr_is_special(x))
    {
        if (fmpr_is_zero(x)) return -1;
        if (fmpr_is_inf(x)) return 1;
        if (fmpr_is_nan(x)) return 0;
        return -1;
    }

    if (fmpz_is_pm1(fmpr_manref(x)))
        return sgn(fmpz_cmp_si(fmpr_expref(x), e));

    bc = fmpz_bits(fmpr_manref(x));

    fmpz_init(t);

    fmpz_add_si_inline(t, fmpr_expref(x), bc);
    fmpz_sub_si_inline(t, t, e);

    ret = (fmpz_sgn(t) <= 0) ? -1 : 1;

    fmpz_clear(t);
    return ret;
}

