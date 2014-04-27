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
arf_mul_special(arf_t z, const arf_t x, const arf_t y)
{
    if (arf_is_zero(x))
    {
        if (arf_is_finite(y))
            arf_zero(z);
        else
            arf_nan(z);
    }
    else if (arf_is_zero(y))
    {
        if (arf_is_finite(x))
            arf_zero(z);
        else
            arf_nan(z);
    }
    else if (arf_is_nan(x) || arf_is_nan(y))
        arf_nan(z);
    else if (arf_sgn(x) == arf_sgn(y))
        arf_pos_inf(z);
    else
        arf_neg_inf(z);
}

