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

/* branch cuts are on (+/-i)*[1,inf] */
int
acb_atan_on_branch_cut(const acb_t z)
{
    arb_t unit;
    int result;

    if (!acb_is_finite(z))
        return 1;

    if (arb_is_nonzero(acb_realref(z)))
        return 0;

    if (arb_contains_si(acb_imagref(z), 1) || arb_contains_si(acb_imagref(z), -1))
        return 1;

    arb_init(unit);
    mag_one(arb_radref(unit));
    result = !arb_contains(unit, acb_imagref(z));
    arb_clear(unit);

    return result;
}

void
acb_atan(acb_t r, const acb_t z, long prec)
{
    if (arb_is_zero(acb_imagref(z)))
    {
        arb_atan(acb_realref(r), acb_realref(z), prec);
        arb_zero(acb_imagref(r));
    }
    else if (!acb_is_finite(z))
    {
        acb_indeterminate(r);
    }
    else
    {
        acb_t t, u;

        acb_init(t);
        acb_init(u);

        if (acb_atan_on_branch_cut(z))
        {
            acb_mul_onei(u, z);
            acb_neg(t, u);
            acb_log1p(t, t, prec);
            acb_log1p(u, u, prec);
            acb_sub(t, t, u, prec);
        }
        else
        {
            acb_onei(t);
            acb_sub(t, t, z, prec);
            acb_div(t, z, t, prec);
            acb_mul_2exp_si(t, t, 1);
            acb_log1p(t, t, prec);
        }

        acb_mul_onei(t, t);
        acb_mul_2exp_si(r, t, -1);

        acb_clear(t);
        acb_clear(u);
    }
}

