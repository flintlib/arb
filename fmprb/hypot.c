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

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include "fmprb.h"

void
fmprb_hypot(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec)
{
    if (fmprb_is_zero(y))
    {
        fmprb_abs(z, x);
    }
    else if (fmprb_is_zero(x))
    {
        fmprb_abs(z, y);
    }
    else
    {
        fmprb_t t;
        fmprb_init(t);
        fmprb_mul(t, x, x, prec + 4);
        fmprb_mul(z, y, y, prec + 4);
        fmprb_add(t, t, z, prec + 4);
        fmprb_sqrtpos(z, t, prec);
        fmprb_clear(t);
    }
}

