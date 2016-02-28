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
_arf_get_integer_mpn(mp_ptr y, mp_srcptr x, mp_size_t xn, slong exp)
{
    slong bot_exp = exp - xn * FLINT_BITS;

    if (bot_exp >= 0)
    {
        mp_size_t bot_limbs;
        mp_bitcnt_t bot_bits;

        bot_limbs = bot_exp / FLINT_BITS;
        bot_bits = bot_exp % FLINT_BITS;

        flint_mpn_zero(y, bot_limbs);

        if (bot_bits == 0)
            flint_mpn_copyi(y + bot_limbs, x, xn);
        else
            y[bot_limbs + xn] = mpn_lshift(y + bot_limbs, x, xn, bot_bits);

        /* exact */
        return 0;
    }
    else if (exp <= 0)
    {
        /* inexact */
        return 1;
    }
    else
    {
        mp_size_t top_limbs;
        mp_bitcnt_t top_bits;
        mp_limb_t cy;

        top_limbs = exp / FLINT_BITS;
        top_bits = exp % FLINT_BITS;

        if (top_bits == 0)
        {
            flint_mpn_copyi(y, x + xn - top_limbs, top_limbs);
            /* inexact */
            return 1;
        }
        else
        {
            /* can be inexact */
            cy = mpn_rshift(y, x + xn - top_limbs - 1,
                top_limbs + 1, FLINT_BITS - top_bits);

            return (cy != 0) || (top_limbs + 1 != xn);
        }
    }
}

