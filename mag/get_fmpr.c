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

#include "mag.h"

void
mag_get_fmpr(fmpr_t x, const mag_t r)
{
    if (mag_is_zero(r))
    {
        fmpr_zero(x);
    }
    else if (mag_is_inf(r))
    {
        fmpr_pos_inf(x);
    }
    else
    {
        fmpr_set_ui_2exp_si(x, MAG_MAN(r), -MAG_BITS);
        _fmpz_add2_fast(fmpr_expref(x), fmpr_expref(x), MAG_EXPREF(r), 0);
    }
}

