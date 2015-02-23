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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "acb.h"

static int
close_to_one(const acb_t z)
{
    mp_limb_t top;

    if (arf_abs_bound_lt_2exp_si(arb_midref(acb_imagref(z))) > -3)
        return 0;

    if (ARF_EXP(arb_midref(acb_realref(z))) == 0)
    {
        ARF_GET_TOP_LIMB(top, arb_midref(acb_realref(z)));

        return (top >> (FLINT_BITS - 4)) == 15;
    }
    else if (ARF_EXP(arb_midref(acb_realref(z))) == 1)
    {
        ARF_GET_TOP_LIMB(top, arb_midref(acb_realref(z)));

        return (top >> (FLINT_BITS - 4)) == 8;
    }

    return 0;
}

void
acb_log(acb_t r, const acb_t z, long prec)
{
#define a acb_realref(z)
#define b acb_imagref(z)

    if (arb_is_zero(b))
    {
        if (arb_is_positive(a))
        {
            arb_log(acb_realref(r), a, prec);
            arb_zero(acb_imagref(r));
        }
        else if (arb_is_negative(a))
        {
            arb_neg(acb_realref(r), a);
            arb_log(acb_realref(r), acb_realref(r), prec);
            arb_const_pi(acb_imagref(r), prec);
        }
        else
        {
            acb_indeterminate(r);
        }
    }
    else if (arb_is_zero(a))
    {
        if (arb_is_positive(b))
        {
            arb_log(acb_realref(r), b, prec);
            arb_const_pi(acb_imagref(r), prec);
            arb_mul_2exp_si(acb_imagref(r), acb_imagref(r), -1);
        }
        else if (arb_is_negative(b))
        {
            arb_neg(acb_realref(r), b);
            arb_log(acb_realref(r), acb_realref(r), prec);
            arb_const_pi(acb_imagref(r), prec);
            arb_mul_2exp_si(acb_imagref(r), acb_imagref(r), -1);
            arb_neg(acb_imagref(r), acb_imagref(r));
        }
        else
        {
            acb_indeterminate(r);
        }
    }
    else
    {
        arb_t t, u;

        arb_init(t);
        arb_init(u);

        if (close_to_one(z))
        {
            arb_sub_ui(u, a, 1, prec + 8);
            arb_mul(t, u, u, prec + 8);
            arb_addmul(t, b, b, prec + 8);
            arb_mul_2exp_si(u, u, 1);
            arb_add(t, t, u, prec + 8);

            arb_log1p(t, t, prec);
            arb_mul_2exp_si(t, t, -1);
        }
        else
        {
            arb_mul(t, a, a, prec + 8);
            arb_addmul(t, b, b, prec + 8);

            if (arb_contains_zero(t) || arf_sgn(arb_midref(t)) < 0)
                arb_zero_pm_inf(t);
            else
                arb_log(t, t, prec);

            arb_mul_2exp_si(t, t, -1);
        }

        acb_arg(u, z, prec);

        arb_swap(acb_realref(r), t);
        arb_swap(acb_imagref(r), u);

        arb_clear(t);
        arb_clear(u);
    }

    if (!acb_is_finite(r))
        acb_indeterminate(r);

#undef a
#undef b
}

