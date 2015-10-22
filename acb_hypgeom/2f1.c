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

#include "acb_hypgeom.h"

static void 
_acb_hypgeom_2f1r_reduced(acb_t res,
    const acb_t b, const acb_t c, const acb_t z, long prec)
{
    acb_t t, u;
    acb_init(t);
    acb_init(u);
    acb_sub_ui(t, z, 1, prec);
    acb_neg(t, t);
    acb_neg(u, b);
    acb_pow(t, t, u, prec);
    acb_rgamma(u, c, prec);
    acb_mul(t, t, u, prec);
    acb_set(res, t);
    acb_clear(t);
    acb_clear(u);
    return;
}

void
acb_hypgeom_2f1(acb_t res, const acb_t a, const acb_t b,
        const acb_t c, const acb_t z, int regularized, long prec)
{
    int algorithm;

    if (!acb_is_finite(a) || !acb_is_finite(b) || !acb_is_finite(c) || !acb_is_finite(z))
    {
        acb_indeterminate(res);
        return;
    }

    if (acb_is_zero(z))
    {
        if (regularized)
            acb_rgamma(res, c, prec);
        else
            acb_one(res);
        return;
    }

    if (regularized && acb_is_int(c) && arb_is_nonpositive(acb_realref(c)))
    {
        if ((acb_is_int(a) && arb_is_nonpositive(acb_realref(a)) &&
            arf_cmp(arb_midref(acb_realref(a)), arb_midref(acb_realref(c))) >= 0) ||
            (acb_is_int(b) && arb_is_nonpositive(acb_realref(b)) &&
            arf_cmp(arb_midref(acb_realref(b)), arb_midref(acb_realref(c))) >= 0))
        {
            acb_zero(res);
            return;
        }
    }

    if (regularized && acb_eq(a, c))
    {
        _acb_hypgeom_2f1r_reduced(res, b, c, z, prec);
        return;
    }

    if (regularized && acb_eq(b, c))
    {
        _acb_hypgeom_2f1r_reduced(res, a, c, z, prec);
        return;
    }

    /* polynomial */
    if (acb_is_int(a) && arf_sgn(arb_midref(acb_realref(a))) <= 0 &&
         arf_cmpabs_ui(arb_midref(acb_realref(a)), prec) < 0)
    {
        acb_hypgeom_2f1_direct(res, a, b, c, z, regularized, prec);
        return;
    }

    /* polynomial */
    if (acb_is_int(b) && arf_sgn(arb_midref(acb_realref(b))) <= 0 &&
         arf_cmpabs_ui(arb_midref(acb_realref(b)), prec) < 0)
    {
        acb_hypgeom_2f1_direct(res, a, b, c, z, regularized, prec);
        return;
    }

    /* Try to reduce to a polynomial case using the Pfaff transformation */
    if (acb_is_exact(c))
    {
        acb_t t;
        acb_init(t);

        acb_sub(t, c, b, prec);

        if (acb_is_int(t) && arb_is_nonpositive(acb_realref(t)))
        {
            acb_hypgeom_2f1_transform(res, a, b, c, z, regularized, 1, prec);
            acb_clear(t);
            return;
        }

        acb_sub(t, c, a, prec);

        if (acb_is_int(t) && arb_is_nonpositive(acb_realref(t)))
        {
            acb_hypgeom_2f1_transform(res, b, a, c, z, regularized, 1, prec);
            acb_clear(t);
            return;
        }

        acb_clear(t);
    }

    /* special value at z = 1 */
    if (acb_is_one(z))
    {
        acb_t t, u, v;

        acb_init(t);
        acb_init(u);
        acb_init(v);

        acb_sub(t, c, a, prec);
        acb_sub(u, c, b, prec);
        acb_sub(v, t, b, prec);

        if (arb_is_positive(acb_realref(v)))
        {
            acb_rgamma(t, t, prec);
            acb_rgamma(u, u, prec);
            acb_mul(t, t, u, prec);
            acb_gamma(v, v, prec);
            acb_mul(t, t, v, prec);

            if (!regularized)
            {
                acb_gamma(v, c, prec);
                acb_mul(t, t, v, prec);
            }

            acb_set(res, t);
        }
        else
        {
            acb_indeterminate(res);
        }

        acb_clear(t);
        acb_clear(u);
        acb_clear(v);

        return;
    }

    algorithm = acb_hypgeom_2f1_choose(z);

    if (algorithm == 0)
    {
        acb_hypgeom_2f1_direct(res, a, b, c, z, regularized, prec);
    }
    else if (algorithm >= 1 && algorithm <= 5)
    {
        acb_hypgeom_2f1_transform(res, a, b, c, z, regularized, algorithm, prec);
    }
    else
    {
        acb_hypgeom_2f1_corner(res, a, b, c, z, regularized, prec);
    }
}

