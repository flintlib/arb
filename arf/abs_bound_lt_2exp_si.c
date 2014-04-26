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

long
arf_abs_bound_lt_2exp_si(const arf_t x)
{
    long res;

    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
            return -ARF_PREC_EXACT;
        else
            return ARF_PREC_EXACT;
    }

    if (!COEFF_IS_MPZ(ARF_EXP(x)))
        return ARF_EXP(x);

    if (fmpz_fits_si(ARF_EXPREF(x)))
        res = fmpz_get_si(ARF_EXPREF(x));
    else
        res = fmpz_sgn(ARF_EXPREF(x)) < 0 ? -ARF_PREC_EXACT : ARF_PREC_EXACT;

    if (res < -ARF_PREC_EXACT)
        res = -ARF_PREC_EXACT;
    if (res > ARF_PREC_EXACT)
        res = ARF_PREC_EXACT;

    return res;
}

