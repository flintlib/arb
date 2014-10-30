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

void
arf_ceil(arf_t z, const arf_t x)
{
    if (arf_is_special(x) || arf_is_int(x))
    {
        arf_set(z, x);
    }
    else
    {
        long exp = ARF_EXP(x);

        /* now exp cannot be too large, as we would have
           caught this in arf_is_int() */
        if (COEFF_IS_MPZ(exp) || exp <= 0)
        {
            if (ARF_SGNBIT(x))
                arf_zero(z);
            else
                arf_one(z);
        }
        else if (exp == 1)
        {
            arf_set_si(z, ARF_SGNBIT(x) ? -1 : 2);
        }
        else
        {
            arf_set_round(z, x, exp, ARF_RND_CEIL);
        }
    }
}

