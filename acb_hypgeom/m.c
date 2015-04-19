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
acb_hypgeom_m_asymp(acb_t res, const acb_t a, const acb_t b, const acb_t z, int regularized, long prec)
{
    acb_t t, u, v, c;

    acb_init(t);
    acb_init(u);
    acb_init(v);
    acb_init(c);

    acb_sub(c, b, a, prec);
    acb_neg(v, z);

    acb_hypgeom_u_asymp(t, a, b, z, -1, prec);
    acb_hypgeom_u_asymp(u, c, b, v, -1, prec);

    /* gamma(b-a) */
    acb_rgamma(v, c, prec);
    acb_mul(t, t, v, prec);

    /* z^(a-b) */
    acb_neg(c, c);
    acb_pow(v, z, c, prec);
    acb_mul(u, u, v, prec);

    /* gamma(a) */
    acb_rgamma(v, a, prec);
    acb_mul(u, u, v, prec);

    /* exp(z) */
    acb_exp(v, z, prec);
    acb_mul(u, u, v, prec);

    /* (-z)^(-a) */
    acb_neg(c, a);
    acb_neg(v, z);
    acb_pow(v, v, c, prec);
    acb_mul(t, t, v, prec);

    acb_add(t, t, u, prec);

    if (!regularized)
    {
        acb_gamma(v, b, prec);
        acb_mul(t, t, v, prec);
    }

    acb_swap(res, t);

    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
    acb_clear(c);
}

static void
_acb_hypgeom_m_1f1(acb_t res, const acb_t a, const acb_t b, const acb_t z, long prec)
{
    if (acb_is_one(a))
    {
        acb_hypgeom_pfq_direct(res, NULL, 0, b, 1, z, -1, prec);
    }
    else
    {
        acb_struct c[3];
        c[0] = *a;
        c[1] = *b;

        acb_init(c + 2);
        acb_one(c + 2);

        acb_hypgeom_pfq_direct(res, c, 1, c + 1, 2, z, -1, prec);

        acb_clear(c + 2);
    }
}

void
acb_hypgeom_m_1f1(acb_t res, const acb_t a, const acb_t b, const acb_t z, int regularized, long prec)
{
    acb_t t;

    if (regularized)
    {
        acb_init(t);
        acb_rgamma(t, b, prec);
    }

    if (arf_sgn(arb_midref(acb_realref(z))) >= 0
        || (acb_is_int(a) && arb_is_nonpositive(acb_realref(a))))
    {
        _acb_hypgeom_m_1f1(res, a, b, z, prec);
    }
    else
    {
        /* Kummer's transformation */
        acb_t u, v;
        acb_init(u);
        acb_init(v);

        acb_sub(u, b, a, prec);
        acb_neg(v, z);

        _acb_hypgeom_m_1f1(u, u, b, v, prec);
        acb_exp(v, z, prec);
        acb_mul(res, u, v, prec);

        acb_clear(u);
        acb_clear(v);
    }

    if (regularized)
    {
        acb_mul(res, res, t, prec);
        acb_clear(t);
    }
}

void
acb_hypgeom_m(acb_t res, const acb_t a, const acb_t b, const acb_t z, int regularized, long prec)
{
    long m = LONG_MAX;
    long n = LONG_MAX;

    if (acb_is_int(a) &&
            arf_cmpabs_2exp_si(arb_midref(acb_realref(a)), 30) < 0)
    {
        m = arf_get_si(arb_midref(acb_realref(a)), ARF_RND_DOWN);
    }

    if (acb_is_int(b) &&
            arf_cmpabs_2exp_si(arb_midref(acb_realref(b)), 30) < 0)
    {
        n = arf_get_si(arb_midref(acb_realref(b)), ARF_RND_DOWN);
    }

    /* terminating */
    if (m <= 0 && m < n && m > -10 * prec && (n > 0 || !regularized))
    {
        acb_hypgeom_m_1f1(res, a, b, z, regularized, prec);
        return;
    }

    /* large */
    if (acb_hypgeom_u_use_asymp(z, prec))
    {
        acb_hypgeom_m_asymp(res, a, b, z, regularized, prec);
        return;
    }

    /* remove singularity */
    if (n <= 0 && n > -10 * prec && regularized)
    {
        acb_t c, d, t, u;

        acb_init(c);
        acb_init(d);
        acb_init(t);
        acb_init(u);

        acb_sub(c, a, b, prec);
        acb_add_ui(c, c, 1, prec);

        acb_neg(d, b);
        acb_add_ui(d, d, 2, prec);

        acb_hypgeom_m_1f1(t, c, d, z, 0, prec);

        acb_pow_ui(u, z, 1 - n, prec);
        acb_mul(t, t, u, prec);

        acb_rising_ui(u, a, 1 - n, prec);
        acb_mul(t, t, u, prec);

        arb_fac_ui(acb_realref(u), 1 - n, prec);
        acb_div_arb(res, t, acb_realref(u), prec);

        acb_clear(c);
        acb_clear(d);
        acb_clear(t);
        acb_clear(u);
    }
    else
    {
        acb_hypgeom_m_1f1(res, a, b, z, regularized, prec);
    }
}

