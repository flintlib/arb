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
acb_root_ui(acb_t res, const acb_t z, ulong n, slong prec)
{
    if (n == 0)
    {
        acb_indeterminate(res);
    }
    else if (n == 1)
    {
        acb_set_round(res, z, prec);
    }
    else if (n == 2)
    {
        acb_sqrt(res, z, prec);
    }
    else if (n == 4)
    {
        acb_sqrt(res, z, prec + 4);
        acb_sqrt(res, res, prec);
    }
    else if (acb_is_real(z) && arb_is_nonnegative(acb_realref(z)))
    {
        arb_root(acb_realref(res), acb_realref(z), n, prec);
        arb_zero(acb_imagref(res));
    }
    else
    {
        acb_log(res, z, prec + 4);
        acb_div_ui(res, res, n, prec + 4);
        acb_exp(res, res, prec);
    }
}

