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

    Copyright (C) 2014 Fredrik Johansson

******************************************************************************/

#include "acb_hypgeom.h"

/* IMAG: erf(z) = 2z/sqrt(pi) * 1F1(1/2, 3/2, -z^2) */
void
acb_hypgeom_erf_1f1a(acb_t res, const acb_t z, long prec)
{
    acb_t a, t, w;
    acb_struct b[2];

    acb_init(a);
    acb_init(b);
    acb_init(b + 1);
    acb_init(t);
    acb_init(w);

    acb_one(a);
    acb_mul_2exp_si(a, a, -1);
    acb_set_ui(b, 3);
    acb_mul_2exp_si(b, b, -1);
    acb_one(b + 1);

    acb_mul(w, z, z, prec);
    acb_neg(w, w);

    acb_hypgeom_pfq_direct(t, a, 1, b, 2, w, -1, prec);

    acb_const_pi(w, prec);
    acb_rsqrt(w, w, prec);
    acb_mul(t, t, w, prec);
    acb_mul(t, t, z, prec);

    acb_mul_2exp_si(res, t, 1);

    acb_clear(a);
    acb_clear(b);
    acb_clear(b + 1);
    acb_clear(t);
    acb_clear(w);
}

/* REAL: erf(x) = 2x/sqrt(pi) * exp(-x^2) 1F1(1, 3/2, x^2) */
void
acb_hypgeom_erf_1f1b(acb_t res, const acb_t z, long prec)
{
    acb_t a, b, t, w;

    acb_init(a);
    acb_init(b);
    acb_init(t);
    acb_init(w);

    acb_set_ui(b, 3);
    acb_mul_2exp_si(b, b, -1);

    acb_mul(w, z, z, prec);

    acb_hypgeom_pfq_direct(t, a, 0, b, 1, w, -1, prec);

    acb_neg(w, w);
    acb_exp(w, w, prec);
    acb_mul(t, t, w, prec);

    acb_const_pi(w, prec);
    acb_rsqrt(w, w, prec);
    acb_mul(t, t, w, prec);
    acb_mul(t, t, z, prec);

    acb_mul_2exp_si(res, t, 1);

    acb_clear(a);
    acb_clear(b);
    acb_clear(t);
    acb_clear (w);
}

void
acb_hypgeom_erf_asymp(acb_t res, const acb_t z, long prec, long prec2)
{
    acb_t a, t, u;

    acb_init(a);
    acb_init(t);
    acb_init(u);

    acb_one(a);
    acb_mul_2exp_si(a, a, -1);
    acb_mul(t, z, z, prec2);

    acb_hypgeom_u_asymp(u, a, a, t, -1, prec2);

    acb_neg(t, t);
    acb_exp(t, t, prec2);
    acb_mul(u, u, t, prec2);

    acb_const_pi(t, prec2);
    acb_sqrt(t, t, prec2);
    acb_mul(t, t, z, prec2);

    acb_div(u, u, t, prec2);

    /* branch cut term: -1 or 1 */
    if (arb_contains_zero(acb_realref(z)))
    {
        arb_zero(acb_imagref(t));
        arf_zero(arb_midref(acb_realref(t)));
        mag_one(arb_radref(acb_realref(t)));
    }
    else
    {
        acb_set_si(t, arf_sgn(arb_midref(acb_realref(z))));
    }

    acb_sub(t, t, u, prec);

    if (arb_is_zero(acb_imagref(z)))
        arb_zero(acb_imagref(t));
    else if (arb_is_zero(acb_realref(z)))
        arb_zero(acb_realref(t));

    acb_set(res, t);

    acb_clear(a);
    acb_clear(t);
    acb_clear(u);
}

void
acb_hypgeom_erf(acb_t res, const acb_t z, long prec)
{
    double x, y, absz2, logz;
    long prec2;

    if (!acb_is_finite(z))
    {
        acb_indeterminate(res);
        return;
    }

    if (acb_is_zero(z))
    {
        acb_zero(res);
        return;
    }

    if ((arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), 0) < 0 &&
         arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), 0) < 0))
    {
        acb_hypgeom_erf_1f1a(res, z, prec);
        return;
    }

    if ((arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), 64) > 0 ||
         arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), 64) > 0))
    {
        acb_hypgeom_erf_asymp(res, z, prec, prec);
        return;
    }

    x = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_DOWN);
    y = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_DOWN);

    absz2 = x * x + y * y;
    logz = 0.5 * log(absz2);

    if (logz - absz2 < -(prec + 8) * 0.69314718055994530942)
    {
        /* If the asymptotic term is small, we can
           compute with reduced precision */
        prec2 = prec + 4 + (y*y - x*x - logz) * 1.4426950408889634074;
        prec2 = FLINT_MAX(8, prec2);
        prec2 = FLINT_MIN(prec2, prec);

        acb_hypgeom_erf_asymp(res, z, prec, prec2);
    }
    else if (arf_cmpabs(arb_midref(acb_imagref(z)), arb_midref(acb_realref(z))) > 0)
    {
        acb_hypgeom_erf_1f1a(res, z, prec);
    }
    else
    {
        acb_hypgeom_erf_1f1b(res, z, prec);
    }
}

