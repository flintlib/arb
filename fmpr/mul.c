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

    Copyright (C) 2012, 2014 Fredrik Johansson

******************************************************************************/

#include "fmpr.h"

long
fmpr_mul(fmpr_t z, const fmpr_t x, const fmpr_t y, long prec, fmpr_rnd_t rnd)
{
    fmpz xv, yv;

    if (fmpr_is_special(x) || fmpr_is_special(y))
    {
        if (fmpr_is_zero(x) && fmpr_is_finite(y))
        {
            fmpr_zero(z);
        }
        else if (fmpr_is_zero(y) && fmpr_is_finite(x))
        {
            fmpr_zero(z);
        }
        else if ((fmpr_is_inf(x) && (fmpr_is_inf(y) || !fmpr_is_special(y))) ||
            (fmpr_is_inf(y) && !fmpr_is_special(x)))
        {
            if (fmpr_sgn(x) == fmpr_sgn(y))
                fmpr_pos_inf(z);
            else
                fmpr_neg_inf(z);
        }
        else
        {
            fmpr_nan(z);
        }

        return FMPR_RESULT_EXACT;
    }

    xv = *fmpr_manref(x);
    yv = *fmpr_manref(y);

    if (!COEFF_IS_MPZ(xv) && !COEFF_IS_MPZ(yv))
    {
        return _fmpr_mul_1x1(z, FLINT_ABS(xv), fmpr_expref(x),
            FLINT_ABS(yv), fmpr_expref(y), (xv ^ yv) < 0, prec, rnd);
    }
    else
    {
        long xn, yn;
        int xsign, ysign;
        mp_limb_t xtmp, ytmp;
        mp_ptr xptr, yptr;

        FMPZ_GET_MPN_READONLY(xsign, xn, xptr, xtmp, xv)
        FMPZ_GET_MPN_READONLY(ysign, yn, yptr, ytmp, yv)

        if (xn >= yn)
            return _fmpr_mul_mpn(z, xptr, xn, fmpr_expref(x),
                yptr, yn, fmpr_expref(y), xsign ^ ysign, prec, rnd);
        else
            return _fmpr_mul_mpn(z, yptr, yn, fmpr_expref(y),
                xptr, xn, fmpr_expref(x), ysign ^ xsign, prec, rnd);
    }
}

