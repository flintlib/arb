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

    Copyright (C) 2016 Arb authors

******************************************************************************/

#include "acb_hypgeom.h"

void
acb_hypgeom_gamma_lower_1f1a(acb_t res, const acb_t s,
        const acb_t z, int regularized, slong prec)
{
    acb_t a, w;
    acb_struct b[2];

    acb_init(a);
    acb_init(b);
    acb_init(b + 1);
    acb_init(w);

    acb_set(a, s);
    acb_add_ui(b, s, 1, prec);
    acb_one(b + 1);
    acb_neg(w, z);

    /* res = 1F1(s, s+1, -z) / s */
    acb_hypgeom_pfq_direct(res, a, 1, b, 2, w, -1, prec);
    acb_div(res, res, s, prec);

    if (regularized == 0)
    {
        acb_pow(a, z, s, prec);
        acb_mul(res, res, a, prec);
    }
    else if (regularized == 1)
    {
        acb_pow(a, z, s, prec);
        acb_mul(res, res, a, prec);
        acb_rgamma(a, s, prec);
        acb_mul(res, res, a, prec);
    }

    acb_clear(a);
    acb_clear(b);
    acb_clear(b + 1);
    acb_clear(w);
}

void
acb_hypgeom_gamma_lower_1f1b(acb_t res, const acb_t s,
        const acb_t z, int regularized, slong prec)
{
    acb_t a, b;

    acb_init(a);
    acb_init(b);

    acb_add_ui(b, s, 1, prec);
    acb_hypgeom_pfq_direct(res, NULL, 0, b, 1, z, -1, prec);
    acb_div(res, res, s, prec);

    acb_neg(a, z);
    acb_exp(a, a, prec);
    acb_mul(res, res, a, prec);

    if (regularized == 0)
    {
        acb_pow(a, z, s, prec);
        acb_mul(res, res, a, prec);
    }
    else if (regularized == 1)
    {
        acb_pow(a, z, s, prec);
        acb_mul(res, res, a, prec);
        acb_rgamma(a, s, prec);
        acb_mul(res, res, a, prec);
    }

    acb_clear(a);
    acb_clear(b);
}

void
acb_hypgeom_gamma_lower(acb_t res, const acb_t s, const acb_t z, int regularized, slong prec)
{
    if (!acb_is_finite(s) || !acb_is_finite(z))
    {
        acb_indeterminate(res);
        return;
    }

    /* todo: handle the case where z is zero */

    /* todo: handle the case where s is an integer */

    if (arf_sgn(arb_midref(acb_realref(z))) > 0)
        acb_hypgeom_gamma_lower_1f1b(res, s, z, regularized, prec);
    else
        acb_hypgeom_gamma_lower_1f1a(res, s, z, regularized, prec);
}
