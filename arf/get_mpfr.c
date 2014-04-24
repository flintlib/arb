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
arf_get_mpfr(mpfr_t x, const arf_t y, mpfr_rnd_t rnd)
{
    int r;

    if (arf_is_special(y))
    {
        if (arf_is_zero(y))
            mpfr_set_zero(x, 0);
        else if (arf_is_pos_inf(y))
            mpfr_set_inf(x, 1);
        else if (arf_is_neg_inf(y))
            mpfr_set_inf(x, -1);
        else
            mpfr_set_nan(x);
        r = 0;
    }
    else if (COEFF_IS_MPZ(*ARF_EXPREF(y)))
    {
        printf("exception: exponent too large to convert to mpfr");
        abort();
    }
    else
    {
        __mpfr_struct t;
        mp_size_t n;
        mp_srcptr d;

        ARF_GET_MPN_READONLY(d, n, y);

        t._mpfr_d = (mp_ptr) d;
        t._mpfr_exp = ARF_EXP(y);
        t._mpfr_prec = n * FLINT_BITS;
        t._mpfr_sign = ARF_SGNBIT(y) ? -1 : 1;

        r = mpfr_set(x, &t, rnd);
    }

    return r;
}

