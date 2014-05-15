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

#include "acb.h"

void
acb_tan_lower_halfplane(acb_t r, const acb_t z, long prec, int pi, int cot)
{
    arb_t s, c, t, u, v;
    long wp;

#define a acb_realref(z)
#define b acb_imagref(z)

    arb_init(s);
    arb_init(c);
    arb_init(t);
    arb_init(u);
    arb_init(v);

    wp = prec + 6;

    arb_mul_2exp_si(t, a, 1);
    if (pi)
        arb_sin_cos_pi(s, c, t, wp);
    else
        arb_sin_cos(s, c, t, wp);

    /* t = exp(4b) - 1 */
    if (pi)
    {
        arb_const_pi(t, wp);
        arb_mul(t, t, b, wp);
        arb_mul_2exp_si(t, t, 2);
    }
    else
    {
        arb_mul_2exp_si(t, b, 2);
    }

    arb_expm1(t, t, wp);

    /* u = 2exp(2b) (sqrt would be inaccurate when b is very negative) */
    if (pi)
    {
        arb_const_pi(u, wp);
        arb_mul(u, u, b, wp);
        arb_mul_2exp_si(u, u, 1);
    }
    else
    {
        arb_mul_2exp_si(u, b, 1);
    }
    arb_exp(u, u, wp);
    arb_mul_2exp_si(u, u, 1);

    /* im = (exp(4b) - 1) / (2 cos(2a) exp(2b) + (exp(4b) - 1) + 2) */
    arb_mul(v, c, u, wp);
    if (cot)
        arb_neg(v, v);
    arb_add(v, v, t, wp);
    arb_add_ui(v, v, 2, wp);
    arb_div(acb_imagref(r), t, v, prec);
    if (cot)
        arb_neg(acb_imagref(r), acb_imagref(r));

    /* re = 2 exp(2b) sin(2a) / (...) */
    arb_mul(s, s, u, wp);
    arb_div(acb_realref(r), s, v, prec);

    arb_clear(s);
    arb_clear(c);
    arb_clear(t);
    arb_clear(u);
    arb_clear(v);

#undef a
#undef b
}

void
acb_tan_near_real(acb_t r, const acb_t z, long prec, int pi, int cot)
{
#define a acb_realref(z)
#define b acb_imagref(z)

    arb_t sa, ca, sb, cb;
    long wp;

    arb_init(sa);
    arb_init(ca);
    arb_init(sb);
    arb_init(cb);

    wp = prec + 6;

    if (pi)
    {
        arb_mul_2exp_si(sa, a, 1);
        arb_sin_cos_pi(sa, ca, sa, wp);

        arb_const_pi(sb, wp);
        arb_mul(sb, sb, b, wp);
        arb_mul_2exp_si(sb, sb, 1);
        arb_sinh_cosh(sb, cb, sb, wp);
    }
    else
    {
        arb_mul_2exp_si(sa, a, 1);
        arb_sin_cos(sa, ca, sa, wp);
        arb_mul_2exp_si(sb, b, 1);
        arb_sinh_cosh(sb, cb, sb, wp);
    }

    if (cot)
    {
        arb_sub(ca, ca, cb, wp);
    }
    else
    {
        arb_add(ca, ca, cb, wp);
    }

    arb_div(acb_realref(r), sa, ca, prec);
    arb_div(acb_imagref(r), sb, ca, prec);

    if (cot)
        arb_neg(acb_realref(r), acb_realref(r));

    arb_clear(sa);
    arb_clear(ca);
    arb_clear(sb);
    arb_clear(cb);

#undef a
#undef b
}

void
acb_tan(acb_t r, const acb_t z, long prec)
{
    if (arb_is_zero(acb_imagref(z)))
    {
        arb_tan(acb_realref(r), acb_realref(z), prec);
        arb_zero(acb_imagref(r));
    }
    else if (arb_is_zero(acb_realref(z)))
    {
        arb_tanh(acb_imagref(r), acb_imagref(z), prec);
        arb_zero(acb_realref(r));
    }
    else
    {
        if (arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), 1) < 0)
        {
            acb_tan_near_real(r, z, prec, 0, 0);
        }
        else if (arf_sgn(arb_midref(acb_imagref(z))) < 0)
        {
            acb_tan_lower_halfplane(r, z, prec, 0, 0);
        }
        else
        {
            acb_neg(r, z);
            acb_tan_lower_halfplane(r, r, prec, 0, 0);
            acb_neg(r, r);
        }
    }
}

void
acb_tan_pi(acb_t r, const acb_t z, long prec)
{
    if (arb_is_zero(acb_imagref(z)))
    {
        arb_tan_pi(acb_realref(r), acb_realref(z), prec);
        arb_zero(acb_imagref(r));
    }
    else if (arb_is_zero(acb_realref(z)))
    {
        arb_t t;
        arb_init(t);
        arb_const_pi(t, prec + 4);
        arb_mul(t, acb_imagref(z), t, prec + 4);
        arb_tanh(acb_imagref(r), t, prec);
        arb_zero(acb_realref(r));
        arb_clear(t);
    }
    else
    {
        if (arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), 1) < 0)
        {
            acb_tan_near_real(r, z, prec, 1, 0);
        }
        else if (arf_sgn(arb_midref(acb_imagref(z))) < 0)
        {
            acb_tan_lower_halfplane(r, z, prec, 1, 0);
        }
        else
        {
            acb_neg(r, z);
            acb_tan_lower_halfplane(r, r, prec, 1, 0);
            acb_neg(r, r);
        }
    }
}

