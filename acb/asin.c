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
acb_asin(acb_t res, const acb_t z, long prec)
{
    if (arb_is_zero(acb_realref(z)))
    {
        arb_asinh(acb_imagref(res), acb_imagref(z), prec);
        arb_zero(acb_realref(res));
    }
    else
    {
        acb_t t;
        acb_init(t);

        acb_mul(t, z, z, prec);
        acb_sub_ui(t, t, 1, prec);
        acb_neg(t, t);
        acb_sqrt(t, t, prec);

        if (acb_is_real(z) && acb_is_real(t))
        {
            arb_atan2(acb_realref(res), acb_realref(z), acb_realref(t), prec);
            arb_zero(acb_imagref(res));
        }
        else
        {
            acb_mul_onei(res, z);
            acb_add(res, res, t, prec);
            acb_log(res, res, prec);
            acb_div_onei(res, res);
        }

        acb_clear(t);
    }
}

