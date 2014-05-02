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
arf_sub_special(arf_t z, const arf_t x, const arf_t y, long prec, arf_rnd_t rnd)
{
    if (arf_is_zero(x))
    {
        if (arf_is_zero(y))
        {
            arf_zero(z);
            return 0;
        }
        else
            return arf_neg_round(z, y, prec, rnd);
    }
    else if (arf_is_zero(y))
    {
        return arf_set_round(z, x, prec, rnd);
    }
    else if (arf_is_nan(x) || arf_is_nan(y)
        || (arf_is_pos_inf(x) && arf_is_pos_inf(y))
        || (arf_is_neg_inf(x) && arf_is_neg_inf(y)))
    {
        arf_nan(z);
        return 0;
    }
    else if (arf_is_special(x))
    {
        arf_set(z, x);
        return 0;
    }
    else
    {
        arf_neg(z, y);
        return 0;
    }
}

int
arf_sub(arf_ptr z, arf_srcptr x, arf_srcptr y, long prec, arf_rnd_t rnd)
{
    mp_size_t xn, yn;
    mp_srcptr xptr, yptr;
    long shift;

    if (arf_is_special(x) || arf_is_special(y))
    {
        return arf_sub_special(z, x, y, prec, rnd);
    }

    shift = _fmpz_sub_small(ARF_EXPREF(x), ARF_EXPREF(y));

    ARF_GET_MPN_READONLY(xptr, xn, x);
    ARF_GET_MPN_READONLY(yptr, yn, y);

    if (shift >= 0)
        return _arf_add_mpn(z, xptr, xn, ARF_SGNBIT(x), ARF_EXPREF(x),
                           yptr, yn, ARF_SGNBIT(y) ^ 1, shift, prec, rnd);
    else
        return _arf_add_mpn(z, yptr, yn, ARF_SGNBIT(y) ^ 1, ARF_EXPREF(y),
                           xptr, xn, ARF_SGNBIT(x), -shift, prec, rnd);
}

