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

long
fmpr_abs_bound_lt_2exp_si(const fmpr_t x)
{
    if (fmpr_is_special(x))
    {
        if (fmpr_is_zero(x))
            return -FMPR_PREC_EXACT;
        else
            return FMPR_PREC_EXACT;
    }
    else
    {
        long res;
        fmpz_t t;
        fmpz_init(t);
        fmpr_abs_bound_lt_2exp_fmpz(t, x);

        if (fmpz_fits_si(t))
            res = fmpz_get_si(t);
        else
            res = fmpz_sgn(t) < 0 ? -FMPR_PREC_EXACT : FMPR_PREC_EXACT;

        fmpz_clear(t);

        if (res < -FMPR_PREC_EXACT)
            res = -FMPR_PREC_EXACT;
        if (res > FMPR_PREC_EXACT)
            res = FMPR_PREC_EXACT;

        return res;
    }
}

