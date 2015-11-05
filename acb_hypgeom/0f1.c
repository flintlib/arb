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
acb_hypgeom_0f1_asymp(acb_t res, const acb_t a, const acb_t z, int regularized, long prec)
{
    acb_t t, u, v;
    int neg;

    acb_init(t);
    acb_init(u);
    acb_init(v);

    /* both expansions are correct, but we want the one that works better
       on the real line */
    neg = arf_sgn(arb_midref(acb_realref(z))) < 0;

    if (neg)
        acb_neg(t, z);
    else
        acb_set(t, z);

    acb_sqrt(t, t, prec);
    acb_mul_2exp_si(v, t, 1);
    acb_sub_ui(u, a, 1, prec);

    if (neg)
        acb_hypgeom_bessel_j_asymp(v, u, v, prec);
    else
        acb_hypgeom_bessel_i_asymp(v, u, v, prec);

    acb_neg(u, u);
    acb_pow(t, t, u, prec);
    acb_mul(v, v, t, prec);

    if (!regularized)
    {
        acb_gamma(t, a, prec);
        acb_mul(v, v, t, prec);
    }

    acb_set(res, v);

    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
}

void
acb_hypgeom_0f1_direct(acb_t res, const acb_t a, const acb_t z, int regularized, long prec)
{
    if (regularized)
    {
        if (acb_is_int(a) && arf_sgn(arb_midref(acb_realref(a))) <= 0)
        {
            acb_t t, u;
            acb_init(t);
            acb_init(u);
            acb_sub_ui(t, a, 2, prec);
            acb_neg(t, t);
            acb_sub_ui(u, a, 1, prec);
            acb_neg(u, u);
            acb_pow(u, z, u, prec);
            /* this cannot recurse infinitely, because t will
               either be an exact positive integer, or inexact */
            acb_hypgeom_0f1_direct(res, t, z, regularized, prec);
            acb_mul(res, res, u, prec);
            acb_clear(t);
            acb_clear(u);
        }
        else  /* todo: could skip when a=1 or a=2 */
        {
            acb_t t;
            acb_init(t);
            acb_rgamma(t, a, prec);
            acb_hypgeom_0f1_direct(res, a, z, 0, prec);
            acb_mul(res, res, t, prec);
            acb_clear(t);
        }
    }
    else
    {
        acb_struct bb[2];
        bb[0] = *a;
        acb_init(bb + 1);
        acb_one(bb + 1);
        acb_hypgeom_pfq_direct(res, NULL, 0, bb, 2, z, -1, prec);
        acb_clear(bb + 1);
    }
}

int
acb_hypgeom_0f1_use_asymp(const acb_t z, long prec)
{
    double x, y, c;

    if ((arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), 0) < 0 &&
         arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), 0) < 0))
    {
        return 0;
    }

    if ((arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), 0) > 128 ||
         arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), 0) > 128))
    {
        return 1;
    }

    x = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_DOWN);
    y = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_DOWN);

    c = prec * 0.69314718055994530942;
    c = c * c;
    c = c * c;

    return x * x + y * y > c;
}

void
acb_hypgeom_0f1(acb_t res, const acb_t a, const acb_t z, int regularized, long prec)
{
    if (acb_hypgeom_0f1_use_asymp(z, prec))
        acb_hypgeom_0f1_asymp(res, a, z, regularized, prec);
    else
        acb_hypgeom_0f1_direct(res, a, z, regularized, prec);
}

