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
arf_set_fmpr(arf_t y, const fmpr_t x)
{
    if (fmpr_is_special(x))
    {
        if (fmpr_is_zero(x))
            arf_zero(y);
        else if (fmpr_is_pos_inf(x))
            arf_pos_inf(y);
        else if (fmpr_is_neg_inf(x))
            arf_neg_inf(y);
        else
            arf_nan(y);
    }
    else
    {
        arf_set_fmpz(y, fmpr_manref(x));
        fmpz_add_inline(ARF_EXPREF(y), ARF_EXPREF(y), fmpr_expref(x));
    }
}

