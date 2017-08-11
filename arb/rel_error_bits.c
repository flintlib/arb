/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

slong
arb_rel_error_bits(const arb_t x)
{
    fmpz_t t;
    slong result;

    /* fast path for small exponents */
    if (ARB_IS_LAGOM(x))
    {
        if (mag_is_zero(arb_radref(x)))
            return -ARF_PREC_EXACT;
        else if (arf_is_special(arb_midref(x)))
            return ARF_PREC_EXACT;
        else
            return MAG_EXP(arb_radref(x)) + 1 - ARF_EXP(arb_midref(x));
    }

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

