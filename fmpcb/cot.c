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

#include "fmpcb.h"

void fmpcb_tan_lower_halfplane(fmpcb_t r, const fmpcb_t z, long prec, int pi, int cot);

void fmpcb_tan_near_real(fmpcb_t r, const fmpcb_t z, long prec, int pi, int cot);

void
fmpcb_cot(fmpcb_t r, const fmpcb_t z, long prec)
{
    if (fmprb_is_zero(fmpcb_imagref(z)))
    {
        fmprb_cot(fmpcb_realref(r), fmpcb_realref(z), prec);
        fmprb_zero(fmpcb_imagref(r));
    }
    else if (fmprb_is_zero(fmpcb_realref(z)))
    {
        fmprb_coth(fmpcb_imagref(r), fmpcb_imagref(z), prec);
        fmprb_neg(fmpcb_imagref(r), fmpcb_imagref(r));
        fmprb_zero(fmpcb_realref(r));
    }
    else
    {
        if (fmpr_cmpabs_2exp_si(fmprb_midref(fmpcb_imagref(z)), 1) < 0)
        {
            if (fmpr_cmpabs_2exp_si(fmprb_midref(fmpcb_realref(z)), 1) < 0)
            {
                /* cos(...) - cosh(...) becomes inaccurate near 0 */
                fmpcb_tan(r, z, prec + 4);
                fmpcb_inv(r, r, prec);
            }
            else
            {
                fmpcb_tan_near_real(r, z, prec, 0, 1);
            }
        }
        else if (fmpr_sgn(fmprb_midref(fmpcb_imagref(z))) < 0)
        {
            fmpcb_tan_lower_halfplane(r, z, prec, 0, 1);
        }
        else
        {
            fmpcb_neg(r, z);
            fmpcb_tan_lower_halfplane(r, r, prec, 0, 1);
            fmpcb_neg(r, r);
        }
    }
}

void
fmpcb_cot_pi(fmpcb_t r, const fmpcb_t z, long prec)
{
    if (fmprb_is_zero(fmpcb_imagref(z)))
    {
        fmprb_cot_pi(fmpcb_realref(r), fmpcb_realref(z), prec);
        fmprb_zero(fmpcb_imagref(r));
    }
    else if (fmprb_is_zero(fmpcb_realref(z)))
    {
        fmprb_t t;
        fmprb_init(t);
        fmprb_const_pi(t, prec + 4);
        fmprb_mul(t, fmpcb_imagref(z), t, prec + 4);
        fmprb_coth(fmpcb_imagref(r), t, prec);
        fmprb_neg(fmpcb_imagref(r), fmpcb_imagref(r));
        fmprb_zero(fmpcb_realref(r));
        fmprb_clear(t);
    }
    else
    {
        if (fmpr_cmpabs_2exp_si(fmprb_midref(fmpcb_imagref(z)), 1) < 0)
        {
            if (fmpr_cmpabs_2exp_si(fmprb_midref(fmpcb_realref(z)), 1) < 0)
            {
                /* cos(...) - cosh(...) becomes inaccurate near 0 */
                fmpcb_tan_pi(r, z, prec + 4);
                fmpcb_inv(r, r, prec);
            }
            else
            {
                fmpcb_tan_near_real(r, z, prec, 1, 1);
            }
        }
        else if (fmpr_sgn(fmprb_midref(fmpcb_imagref(z))) < 0)
        {
            fmpcb_tan_lower_halfplane(r, z, prec, 1, 1);
        }
        else
        {
            fmpcb_neg(r, z);
            fmpcb_tan_lower_halfplane(r, r, prec, 1, 1);
            fmpcb_neg(r, r);
        }
    }
}

