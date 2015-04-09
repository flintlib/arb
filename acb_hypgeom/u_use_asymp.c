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

#include "acb_hypgeom.h"

int
acb_hypgeom_u_use_asymp(const acb_t z, long prec)
{
    double x, y;

    if ((arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), 0) < 0 &&
         arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), 0) < 0))
    {
        return 0;
    }

    if ((arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), 0) > 64 ||
         arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), 0) > 64))
    {
        return 1;
    }

    x = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_DOWN);
    y = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_DOWN);

    return sqrt(x * x + y * y) > prec * 0.69314718055994530942;
}

