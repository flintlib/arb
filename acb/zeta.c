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

#include "acb.h"
#include "acb_poly.h"

void
acb_zeta_si(acb_t z, slong s, slong prec)
{
    if (s >= 0)
    {
        arb_zeta_ui(acb_realref(z), s, prec);
    }
    else
    {
        arb_bernoulli_ui(acb_realref(z), 1-s, prec);
        arb_div_ui(acb_realref(z), acb_realref(z), 1-s, prec);
        arb_neg(acb_realref(z), acb_realref(z));
    }

    arb_zero(acb_imagref(z));
    return;
}

void
acb_hurwitz_zeta(acb_t z, const acb_t s, const acb_t a, slong prec)
{
    if (acb_is_one(a) && acb_is_int(s) &&
        arf_cmpabs_2exp_si(arb_midref(acb_realref(s)), FLINT_BITS - 1) < 0)
    {
        acb_zeta_si(z, arf_get_si(arb_midref(acb_realref(s)), ARF_RND_DOWN), prec);
        return;
    }

    _acb_poly_zeta_cpx_series(z, s, a, 0, 1, prec);
}

void
acb_zeta(acb_t z, const acb_t s, slong prec)
{
    acb_t a;
    acb_init(a);
    acb_one(a);

    if (acb_is_int(s) &&
        arf_cmpabs_2exp_si(arb_midref(acb_realref(s)), FLINT_BITS - 1) < 0)
    {
        acb_zeta_si(z, arf_get_si(arb_midref(acb_realref(s)), ARF_RND_DOWN), prec);
        return;
    }

    if (arf_sgn(arb_midref(acb_realref(s))) < 0)
    {
        acb_t t, u, v;
        slong wp = prec + 6;

        acb_init(t);
        acb_init(u);
        acb_init(v);

        acb_sub_ui(t, s, 1, wp);

        /* 2 * (2pi)^(s-1) */
        arb_const_pi(acb_realref(u), wp);
        acb_mul_2exp_si(u, u, 1);
        acb_pow(u, u, t, wp);
        acb_mul_2exp_si(u, u, 1);

        /* sin(pi*s/2) */
        acb_mul_2exp_si(v, s, -1);
        acb_sin_pi(v, v, wp);
        acb_mul(u, u, v, wp);

        /* gamma(1-s) zeta(1-s) */
        acb_neg(t, t);
        acb_gamma(v, t, wp);
        acb_mul(u, u, v, wp);
        acb_hurwitz_zeta(v, t, a, wp);
        acb_mul(z, u, v, prec);

        acb_clear(t);
        acb_clear(u);
        acb_clear(v);
    }
    else
    {
        acb_hurwitz_zeta(z, s, a, prec);
    }

    acb_clear(a);
}

