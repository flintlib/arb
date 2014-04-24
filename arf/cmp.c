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
arf_cmp(const arf_t x, const arf_t y)
{
    int xs, ys, ec, mc;
    mp_size_t xn, yn;
    mp_srcptr xp, yp;

    if (arf_is_special(x) || arf_is_special(y))
    {
        if (arf_equal(x, y))
            return 0;
        if (arf_is_nan(x) || arf_is_nan(y))
            return 0;
        if (arf_is_zero(y)) return arf_sgn(x);
        if (arf_is_zero(x)) return -arf_sgn(y);
        if (arf_is_pos_inf(x)) return 1;
        if (arf_is_neg_inf(y)) return 1;
        return -1;
    }

    xs = ARF_SGNBIT(x);
    ys = ARF_SGNBIT(y);

    if (xs != ys)
        return xs ? -1 : 1;

    ec = fmpz_cmp(ARF_EXPREF(x), ARF_EXPREF(y));

    if (ec != 0)
        return ((ec < 0) ^ xs) ? -1 : 1;

    ARF_GET_MPN_READONLY(xp, xn, x);
    ARF_GET_MPN_READONLY(yp, yn, y);

    if (xn >= yn)
        mc = mpn_cmp(xp + xn - yn, yp, yn);
    else
        mc = mpn_cmp(xp, yp + yn - xn, xn);

    if (mc != 0)
        return ((mc < 0) ^ xs) ? -1 : 1;

    if (xn != yn)
        return ((xn < yn) ^ xs) ? -1 : 1;

    return 0;
}

