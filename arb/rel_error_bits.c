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

#include "arb.h"

long
arb_rel_error_bits(const arb_t x)
{
    fmpz_t t;
    long result;

    if (mag_is_zero(arb_radref(x)))
    {
        if (arf_is_nan(arb_midref(x)))
            return ARF_PREC_EXACT;
        else
            return -ARF_PREC_EXACT;
    }

    if (arf_is_special(arb_midref(x)) || mag_is_inf(arb_radref(x)))
        return ARF_PREC_EXACT;

    fmpz_init(t);
    fmpz_add_ui(t, MAG_EXPREF(arb_radref(x)), 1);
    result = _fmpz_sub_small(t, ARF_EXPREF(arb_midref(x)));
    fmpz_clear(t);

    return result;
}

