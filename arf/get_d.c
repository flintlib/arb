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
#include "double_extras.h"

double
arf_get_d(const arf_t x, arf_rnd_t rnd)
{
    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
            return 0.0;
        else if (arf_is_pos_inf(x))
            return D_INF;
        else if (arf_is_neg_inf(x))
            return -D_INF;
        else
            return D_NAN;
    }
    else
    {
        arf_t t;
        mp_srcptr tp;
        mp_size_t tn;
        double v;

        /* also catches bignum exponents */
        if (ARF_EXP(x) > 2000 || ARF_EXP(x) < -2000)
        {
            if (fmpz_sgn(ARF_EXPREF(x)) > 0)
                return ARF_SGNBIT(x) ? -D_INF : D_INF;
            else
                return 0.0;
        }

        arf_init(t);
        arf_set_round(t, x, 53, rnd);
        ARF_GET_MPN_READONLY(tp, tn, t);

        if (tn == 1)
            v = (double)(tp[0]);
        else
            v = (double)(tp[1]) + (double)(tp[0]) * ldexp(1,-32);

        v = ldexp(v, ARF_EXP(t) - FLINT_BITS);

        if (ARF_SGNBIT(t))
            v = -v;

        arf_clear(t);

        return v;
    }
}

