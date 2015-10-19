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

