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

    Copyright (C) 2012, 2013 Fredrik Johansson

******************************************************************************/

#include "arb.h"

void
arb_atan2(arb_t r, const arb_t b, const arb_t a, long prec)
{
#define am arb_midref(a)
#define ar arb_radref(a)
#define bm arb_midref(b)
#define br arb_radref(b)

    /* a + bi is a real number */
    if (arb_is_zero(b))
    {
        /* exactly zero */
        if (arb_is_zero(a))
        {
            /* define arg(0) = 0 by convention */
            arb_zero(r);
        }
        /* interval contains only nonnegative numbers */
        else if (arf_sgn(am) > 0 && arf_cmpabs_mag(am, ar) >= 0)
        {
            arb_zero(r);
        }
        /* interval contains only negative numbers */
        else if (arf_sgn(am) < 0 && arf_cmpabs_mag(am, ar) > 0)
        {
            arb_const_pi(r, prec);
        }
        else
        {
            /* both positive and negative -- argument will be in [0, pi] */
            arb_t t;
            arb_init(t);
            arb_const_pi(t, prec);
            arb_mul_2exp_si(t, t, -1);
            arb_set(r, t);
            arb_add_error(r, t);
            arb_clear(t);
        }
    }
    /* an imaginary number */
    else if (arb_is_zero(a))
    {
        /* interval contains only positive numbers */
        if (arf_sgn(bm) > 0 && arf_cmpabs_mag(bm, br) > 0)
        {
            arb_const_pi(r, prec);
            arb_mul_2exp_si(r, r, -1);
        }
        /* interval contains only negative numbers */
        else if (arf_sgn(bm) < 0 && arf_cmpabs_mag(bm, br) > 0)
        {
            arb_const_pi(r, prec);
            arb_neg(r, r);
            arb_mul_2exp_si(r, r, -1);
        }
        else
        {
            /* both positive and negative -- argument will be in 0 +/- pi/2 */
            arb_t t;
            arb_init(t);
            arb_const_pi(t, MAG_BITS);
            arb_mul_2exp_si(t, t, -1);
            arb_zero(r);
            arb_add_error(r, t);
            arb_clear(t);
        }
    }
    /* strictly in the right half-plane -- atan(b/a) */
    else if (arf_sgn(am) > 0 && arf_cmpabs_mag(am, ar) > 0)
    {
        arb_div(r, b, a, prec);
        arb_atan(r, r, prec);
    }
    /* strictly in the upper half-plane -- pi/2 - atan(a/b) */
    else if (arf_sgn(bm) > 0 && arf_cmpabs_mag(bm, br) > 0)
    {
        arb_t t;
        arb_init(t);
        arb_div(r, a, b, prec);
        arb_atan(r, r, prec);
        arb_const_pi(t, prec);
        arb_mul_2exp_si(t, t, -1);
        arb_sub(r, t, r, prec);
        arb_clear(t);
    }
    /* strictly in the lower half-plane -- -pi/2 - atan(a/b) */
    else if (arf_sgn(bm) < 0 && arf_cmpabs_mag(bm, br) > 0)
    {
        arb_t t;
        arb_init(t);
        arb_div(r, a, b, prec);
        arb_atan(r, r, prec);
        arb_const_pi(t, prec);
        arb_mul_2exp_si(t, t, -1);
        arb_add(r, t, r, prec);
        arb_neg(r, r);
        arb_clear(t);
    }
    /* overlaps the nonpositive half-axis -- [-pi, pi] */
    else
    {
        arb_t t;
        arb_init(t);
        arb_const_pi(t, MAG_BITS);
        arb_zero(r);
        arb_add_error(r, t);
        arb_clear(t);
    }
}

