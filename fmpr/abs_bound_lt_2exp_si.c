/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpr.h"

slong
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
        slong res;
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

