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

#include "fmpcb.h"

void
fmpcb_arg(fmprb_t r, const fmpcb_t z, long prec)
{
#define a fmpcb_realref(z)
#define b fmpcb_imagref(z)
#define am fmprb_midref(a)
#define ar fmprb_radref(a)
#define bm fmprb_midref(b)
#define br fmprb_radref(b)

    /* a real number */
    if (fmprb_is_zero(b))
    {
        /* exactly zero */
        if (fmprb_is_zero(a))
        {
            /* define arg(0) = 0 by convention */
            fmprb_zero(r);
        }
        /* interval contains only nonnegative numbers */
        else if (fmpr_sgn(am) > 0 && fmpr_cmpabs(am, ar) >= 0)
        {
            fmprb_zero(r);
        }
        /* interval contains only positive numbers */
        else if (fmpr_sgn(am) < 0 && fmpr_cmpabs(am, ar) > 0)
        {
            fmprb_const_pi(r, prec);
        }
        else
        {
            /* both positive and negative -- argument will be in [0, pi] */
            fmprb_t t;
            fmprb_init(t);
            fmprb_const_pi(t, prec);
            fmprb_mul_2exp_si(t, t, -1);
            fmprb_zero(r);
            fmprb_add_error(r, t);
            fmprb_clear(t);
        }
    }
    /* an imaginary number */
    else if (fmprb_is_zero(a))
    {
        /* interval contains only positive numbers */
        if (fmpr_sgn(bm) > 0 && fmpr_cmpabs(bm, br) > 0)
        {
            fmprb_const_pi(r, prec);
            fmprb_mul_2exp_si(r, r, -1);
        }
        /* interval contains only negative numbers */
        else if (fmpr_sgn(bm) < 0 && fmpr_cmpabs(bm, br) > 0)
        {
            fmprb_const_pi(r, prec);
            fmprb_neg(r, r);
            fmprb_mul_2exp_si(r, r, -1);
        }
        else
        {
            /* both positive and negative -- argument will be in 0 +/- pi/2 */
            fmprb_t t;
            fmprb_init(t);
            fmprb_const_pi(t, FMPRB_RAD_PREC);
            fmprb_mul_2exp_si(t, t, -1);
            fmprb_zero(r);
            fmprb_add_error(r, t);
            fmprb_clear(t);
        }
    }
    /* strictly in the right half-plane -- atan(b/a) */
    else if (fmpr_sgn(am) > 0 && fmpr_cmpabs(am, ar) > 0)
    {
        fmprb_div(r, b, a, prec);
        fmprb_atan(r, r, prec);
    }
    /* strictly in the upper half-plane -- pi/2 - atan(a/b) */
    else if (fmpr_sgn(bm) > 0 && fmpr_cmpabs(bm, br) > 0)
    {
        fmprb_t t;
        fmprb_init(t);
        fmprb_div(r, a, b, prec);
        fmprb_atan(r, r, prec);
        fmprb_const_pi(t, prec);
        fmprb_mul_2exp_si(t, t, -1);
        fmprb_sub(r, t, r, prec);
        fmprb_clear(t);
    }
    /* strictly in the lower half-plane -- -pi/2 - atan(a/b) */
    else if (fmpr_sgn(bm) < 0 && fmpr_cmpabs(bm, br) > 0)
    {
        fmprb_t t;
        fmprb_init(t);
        fmprb_div(r, a, b, prec);
        fmprb_atan(r, r, prec);
        fmprb_const_pi(t, prec);
        fmprb_mul_2exp_si(t, t, -1);
        fmprb_add(r, t, r, prec);
        fmprb_neg(r, r);
        fmprb_clear(t);
    }
    /* overlaps the nonpositive half-axis -- [-pi, pi] */
    else
    {
        fmprb_t t;
        fmprb_init(t);
        fmprb_const_pi(t, FMPRB_RAD_PREC);
        fmprb_zero(r);
        fmprb_add_error(r, t);
        fmprb_clear(t);
    }
}

