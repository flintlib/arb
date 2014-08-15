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
acb_get_mag(mag_t u, const acb_t z)
{
    if (arb_is_zero(acb_imagref(z)))
    {
        arb_get_mag(u, acb_realref(z));
    }
    else if (arb_is_zero(acb_realref(z)))
    {
        arb_get_mag(u, acb_imagref(z));
    }
    else
    {
        mag_t v;
        mag_init(v);

        arb_get_mag(u, acb_realref(z));
        arb_get_mag(v, acb_imagref(z));

        mag_mul(u, u, u);
        mag_addmul(u, v, v);
        mag_sqrt(u, u);

        mag_clear(v);
    }
}

