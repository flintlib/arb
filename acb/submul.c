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

#include "acb.h"

void
acb_submul(acb_t z, const acb_t x, const acb_t y, long prec)
{
    if (arb_is_zero(acb_imagref(y)))
    {
        arb_submul(acb_imagref(z), acb_imagref(x), acb_realref(y), prec);
        arb_submul(acb_realref(z), acb_realref(x), acb_realref(y), prec);
    }
    else if (arb_is_zero(acb_imagref(x)))
    {
        arb_submul(acb_imagref(z), acb_imagref(y), acb_realref(x), prec);
        arb_submul(acb_realref(z), acb_realref(y), acb_realref(x), prec);
    }
    else
    {
        acb_t t;
        acb_init(t);
        acb_mul(t, x, y, prec);
        acb_sub(z, z, t, prec);
        acb_clear(t);
    }
}

