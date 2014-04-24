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
arf_get_fmpr(fmpr_t y, const arf_t x)
{
    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
            fmpr_zero(y);
        else if (arf_is_pos_inf(x))
            fmpr_pos_inf(y);
        else if (arf_is_neg_inf(x))
            fmpr_neg_inf(y);
        else
            fmpr_nan(y);
    }
    else
    {
        mp_srcptr xptr;
        mp_size_t xn;
        long shift;

        ARF_GET_MPN_READONLY(xptr, xn, x);

        _fmpr_set_round_mpn(&shift, fmpr_manref(y),
            xptr, xn, ARF_SGNBIT(x), FMPR_PREC_EXACT, FMPR_RND_DOWN);

        fmpz_add_si(fmpr_expref(y), ARF_EXPREF(x), shift - xn * FLINT_BITS);
    }
}

