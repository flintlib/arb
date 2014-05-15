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

#include "acb.h"

void acb_tan_lower_halfplane(acb_t r, const acb_t z, long prec, int pi, int cot);

void acb_tan_near_real(acb_t r, const acb_t z, long prec, int pi, int cot);

void
acb_cot(acb_t r, const acb_t z, long prec)
{
    if (arb_is_zero(acb_imagref(z)))
    {
        arb_cot(acb_realref(r), acb_realref(z), prec);
        arb_zero(acb_imagref(r));
    }
    else if (arb_is_zero(acb_realref(z)))
    {
        arb_coth(acb_imagref(r), acb_imagref(z), prec);
        arb_neg(acb_imagref(r), acb_imagref(r));
        arb_zero(acb_realref(r));
    }
    else
    {
        if (arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), 1) < 0)
        {
            if (arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), 1) < 0)
            {
                /* cos(...) - cosh(...) becomes inaccurate near 0 */
                acb_tan(r, z, prec + 4);
                acb_inv(r, r, prec);
            }
            else
            {
                acb_tan_near_real(r, z, prec, 0, 1);
            }
        }
        else if (arf_sgn(arb_midref(acb_imagref(z))) < 0)
        {
            acb_tan_lower_halfplane(r, z, prec, 0, 1);
        }
        else
        {
            acb_neg(r, z);
            acb_tan_lower_halfplane(r, r, prec, 0, 1);
            acb_neg(r, r);
        }
    }
}

void
acb_cot_pi(acb_t r, const acb_t z, long prec)
{
    if (arb_is_zero(acb_imagref(z)))
    {
        arb_cot_pi(acb_realref(r), acb_realref(z), prec);
        arb_zero(acb_imagref(r));
    }
    else if (arb_is_zero(acb_realref(z)))
    {
        arb_t t;
        arb_init(t);
        arb_const_pi(t, prec + 4);
        arb_mul(t, acb_imagref(z), t, prec + 4);
        arb_coth(acb_imagref(r), t, prec);
        arb_neg(acb_imagref(r), acb_imagref(r));
        arb_zero(acb_realref(r));
        arb_clear(t);
    }
    else
    {
        if (arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), 1) < 0)
        {
            if (arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), 1) < 0)
            {
                /* cos(...) - cosh(...) becomes inaccurate near 0 */
                acb_tan_pi(r, z, prec + 4);
                acb_inv(r, r, prec);
            }
            else
            {
                acb_tan_near_real(r, z, prec, 1, 1);
            }
        }
        else if (arf_sgn(arb_midref(acb_imagref(z))) < 0)
        {
            acb_tan_lower_halfplane(r, z, prec, 1, 1);
        }
        else
        {
            acb_neg(r, z);
            acb_tan_lower_halfplane(r, r, prec, 1, 1);
            acb_neg(r, r);
        }
    }
}

