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

    Copyright (C) 2015 Fredrik Johansson

******************************************************************************/

#include "acb.h"

void
acb_acosh(acb_t res, const acb_t z, slong prec)
{
    if (acb_is_one(z))
    {
        acb_zero(res);
    }
    else
    {
        acb_t t, u;
        acb_init(t);
        acb_init(u);

        acb_add_ui(t, z, 1, prec);
        acb_sub_ui(u, z, 1, prec);
        acb_sqrt(t, t, prec);
        acb_sqrt(u, u, prec);
        acb_mul(t, t, u, prec);
        acb_add(t, t, z, prec);
        acb_log(res, t, prec);

        acb_clear(t);
        acb_clear(u);
    }
}

