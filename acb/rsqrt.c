/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb.h"

static void
_acb_rsqrt1(acb_t y, const acb_t x, slong prec)
{
#define a acb_realref(x)
#define b acb_imagref(x)
#define c acb_realref(y)
#define d acb_imagref(y)

    arb_t r, t, u, v;
    slong wp;

    /* based on the identity sqrt(z) = sqrt(r) (z+r) / |z+r|: */
    /* 1/sqrt(a+bi) = (1/v)((a+r) - b*i), r = |a+bi|, v = sqrt(r*(b^2+(a+r)^2)) */

    wp = prec + 6;

    arb_init(r);
    arb_init(t);
    arb_init(u);
    arb_init(v);

    /* u = b^2, r = |a+bi| */
    arb_mul(t, a, a, wp);
    arb_mul(u, b, b, wp);
    arb_add(r, t, u, wp);
    arb_sqrtpos(r, r, wp);

    /* t = a+r, v = r*(b^2+(a+r)^2) */
    arb_add(t, r, a, wp);
    arb_mul(v, t, t, wp);
    arb_add(v, v, u, wp);
    arb_mul(v, v, r, wp);

    /* v = 1/sqrt(v) */
    arb_rsqrt(v, v, wp);

    arb_mul(c, t, v, prec);
    arb_mul(d, b, v, prec);
    arb_neg(d, d);

    arb_clear(r);
    arb_clear(t);
    arb_clear(u);
    arb_clear(v);

#undef a
#undef b
#undef c
#undef d
}

static void
_acb_rsqrt(acb_t y, const acb_t x, slong prec)
{
    if (acb_rel_accuracy_bits(x) > 10)
    {
        _acb_rsqrt1(y, x, prec);
    }
    else if (acb_contains_zero(x))
    {
        acb_indeterminate(y);
    }
    else
    {
        /* use derivative (todo: compute better bounds) */
        mag_t t, u;

        mag_init(t);
        mag_init(u);

        acb_get_mag_lower(t, x);

        mag_rsqrt(u, t);
        mag_div(u, u, t);
        mag_hypot(t, arb_radref(acb_realref(x)), arb_radref(acb_imagref(x)));
        mag_mul(u, u, t);

        acb_set(y, x);
        mag_zero(arb_radref(acb_realref(y)));
        mag_zero(arb_radref(acb_imagref(y)));

        _acb_rsqrt1(y, y, prec);
        acb_add_error_mag(y, u);

        mag_clear(t);
        mag_clear(u);
    }
}

void
acb_rsqrt(acb_t y, const acb_t x, slong prec)
{
#define a acb_realref(x)
#define b acb_imagref(x)
#define c acb_realref(y)
#define d acb_imagref(y)

    if (arb_is_zero(b))
    {
        if (arb_is_nonnegative(a))
        {
            arb_rsqrt(c, a, prec);
            arb_zero(d);
            return;
        }
        else if (arb_is_nonpositive(a))
        {
            arb_neg(d, a);
            arb_rsqrt(d, d, prec);
            arb_neg(d, d);
            arb_zero(c);
            return;
        }
    }

    if (arb_is_zero(a))
    {
        if (arb_is_nonnegative(b))
        {
            arb_mul_2exp_si(c, b, 1);
            arb_rsqrt(c, c, prec);
            arb_neg(d, c);
            return;
        }
        else if (arb_is_nonpositive(b))
        {
            arb_mul_2exp_si(c, b, 1);
            arb_neg(c, c);
            arb_rsqrt(c, c, prec);
            arb_set(d, c);
            return;
        }
    }

    if (arb_is_positive(a))
    {
        _acb_rsqrt(y, x, prec);
    }
    else if (arb_is_nonnegative(b))
    {
        acb_neg(y, x);
        _acb_rsqrt(y, y, prec);
        acb_div_onei(y, y);
    }
    else if (arb_is_negative(b))
    {
        acb_neg(y, x);
        _acb_rsqrt(y, y, prec);
        acb_mul_onei(y, y);
    }
    else
    {
        /* todo: more elegant solution? */
        acb_t t;
        acb_init(t);
        acb_neg(y, x);
        _acb_rsqrt(y, y, prec);
        acb_mul_onei(t, y);
        acb_div_onei(y, y);
        acb_union(y, y, t, prec);
        acb_clear(t);
    }

#undef a
#undef b
#undef c
#undef d
}

void
acb_rsqrt_analytic(acb_ptr res, const acb_t z, int analytic, slong prec)
{
    if (analytic && arb_contains_zero(acb_imagref(z)) &&
        !arb_is_positive(acb_realref(z)))
    {
        acb_indeterminate(res);
    }
    else
    {
        acb_rsqrt(res, z, prec);
    }
}

