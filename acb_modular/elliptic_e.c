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

#include "acb_modular.h"

void
acb_modular_elliptic_e(acb_t res, const acb_t m, long prec)
{
    if (acb_is_zero(m))
    {
        acb_const_pi(res, prec);
        acb_mul_2exp_si(res, res, -1);
    }
    else if (acb_is_one(m))
    {
        acb_one(res);
    }
    else
    {
        acb_struct t[2];

        acb_init(t + 0);
        acb_init(t + 1);

        acb_modular_elliptic_k_cpx(t, m, 2, prec);
        acb_mul(t + 1, t + 1, m, prec);
        acb_mul_2exp_si(t + 1, t + 1, 1);
        acb_add(t, t, t + 1, prec);
        acb_sub_ui(t + 1, m, 1, prec);
        acb_mul(res, t, t + 1, prec);
        acb_neg(res, res);

        acb_clear(t + 0);
        acb_clear(t + 1);
    }
}

