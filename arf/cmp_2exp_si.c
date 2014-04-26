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

#include "arf.h"

int
arf_cmp_2exp_si(const arf_t x, long e)
{
    if (arf_is_special(x))
    {
        if (arf_is_zero(x)) return -1;
        if (arf_is_pos_inf(x)) return 1;
        if (arf_is_neg_inf(x)) return -1;
        return 0;
    }

    if (ARF_SGNBIT(x))
        return -1;

    /* Fast path. */
    if (!COEFF_IS_MPZ(ARF_EXP(x)))
    {
        if (ARF_IS_POW2(x) && (ARF_EXP(x) - 1 == e))
            return 0;
        else
            return (ARF_EXP(x) <= e) ? -1 : 1;
    }

    if (ARF_IS_POW2(x))
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_one(t);
        fmpz_add_si(t, t, e);
        if (fmpz_equal(ARF_EXPREF(x), t))
        {
            fmpz_clear(t);
            return 0;
        }
        fmpz_clear(t);
    }

    return (fmpz_cmp_si(ARF_EXPREF(x), e) <= 0) ? -1 : 1;
}

