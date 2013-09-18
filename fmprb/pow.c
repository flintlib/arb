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

    Copyright (C) 2012, 2013 Fredrik Johansson

******************************************************************************/

#include "fmprb.h"

#define BINEXP_LIMIT 64

void
_fmprb_pow_exp(fmprb_t z, const fmprb_t x, int negx, const fmprb_t y, long prec)
{
    fmprb_t t;
    fmprb_init(t);

    if (negx)
    {
        fmprb_neg(t, x);
        fmprb_log(t, t, prec);
    }
    else
        fmprb_log(t, x, prec);

    fmprb_mul(t, t, y, prec);
    fmprb_exp(z, t, prec);
    fmprb_clear(t);
}

void
fmprb_pow(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec)
{
    const fmpr_struct * ymid = fmprb_midref(y);
    const fmpr_struct * yrad = fmprb_radref(y);

    if (fmprb_is_zero(y))
    {
        fmprb_one(z);
        return;
    }

    if (fmpr_is_zero(yrad) && !fmpr_is_special(fmprb_midref(x)))
    {
        /* small half-integer or integer */
        if (fmpr_cmpabs_2exp_si(ymid, BINEXP_LIMIT) < 0 &&
            fmpr_is_int_2exp_si(ymid, -1))
        {
            fmpz_t e;
            fmpz_init(e);            

            if (fmpr_is_int(ymid))
            {
                fmpr_get_fmpz_fixed_si(e, ymid, 0);
                fmprb_pow_fmpz_binexp(z, x, e, prec);
            }
            else
            {
                fmpr_get_fmpz_fixed_si(e, ymid, -1);
                fmprb_sqrt(z, x, prec + fmpz_bits(e));
                fmprb_pow_fmpz_binexp(z, z, e, prec);
            }

            fmpz_clear(e);
            return;
        }
        else if (fmpr_is_int(ymid) && fmpr_sgn(fmprb_midref(x)) < 0)
        {
            /* use (-x)^n = (-1)^n * x^n to avoid NaNs
               at least at high enough precision */
            int odd = !fmpr_is_int_2exp_si(ymid, 1);
            _fmprb_pow_exp(z, x, 1, y, prec);
            if (odd)
                fmprb_neg(z, z);
            return;
        }
    }

    _fmprb_pow_exp(z, x, 0, y, prec);
}

