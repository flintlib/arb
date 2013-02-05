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
fmpr_get_mpfr(mpfr_t x, const fmpr_t y, mpfr_rnd_t rnd)
{
    int r;

    if (fmpr_is_special(y))
    {
        if (fmpr_is_zero(y)) mpfr_set_zero(x, 0);
        else if (fmpr_is_pos_inf(y)) mpfr_set_inf(x, 1);
        else if (fmpr_is_neg_inf(y)) mpfr_set_inf(x, -1);
        else mpfr_set_nan(x);
        r = 0;
    }
    else if (COEFF_IS_MPZ(*fmpr_expref(y)))
    {
        printf("exception: exponent too large to convert to mpfr");
        abort();
    }
    else
    {
        if (!COEFF_IS_MPZ(*fmpr_manref(y)))
            r = mpfr_set_si_2exp(x, *fmpr_manref(y), *fmpr_expref(y), rnd);
        else
            r = mpfr_set_z_2exp(x, COEFF_TO_PTR(*fmpr_manref(y)), *fmpr_expref(y), rnd);

        if (!mpfr_regular_p(x))
        {
            printf("exception: exponent too large to convert to mpfr");
            abort();
        }
    }

    return r;
}

