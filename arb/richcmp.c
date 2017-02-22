/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

int arb_eq(const arb_t x, const arb_t y)
{
    if (arb_is_indeterminate(x) || arb_is_indeterminate(y))
        return 0;

    if (mag_is_inf(arb_radref(x)) || mag_is_inf(arb_radref(y)))
        return 0;

    if (arf_is_inf(arb_midref(x)) || arf_is_inf(arb_midref(y)) ||
        (mag_is_zero(arb_radref(x)) && mag_is_zero(arb_radref(y))))
        return arf_equal(arb_midref(x), arb_midref(y));

    return 0;
}

int arb_ne(const arb_t x, const arb_t y)
{
    if (arb_is_indeterminate(x) || arb_is_indeterminate(y))
        return 0;

    if (mag_is_inf(arb_radref(x)) || mag_is_inf(arb_radref(y)))
        return 0;

    if (arf_is_inf(arb_midref(x)) || arf_is_inf(arb_midref(y)) ||
        (mag_is_zero(arb_radref(x)) && mag_is_zero(arb_radref(y))))
        return !arf_equal(arb_midref(x), arb_midref(y));

    return !arb_overlaps(x, y);
}

int arb_lt(const arb_t x, const arb_t y)
{
    arf_struct u[4];
    arf_t t;
    int res;

    if (arb_is_indeterminate(x) || arb_is_indeterminate(y))
        return 0;

    if (mag_is_inf(arb_radref(x)) || mag_is_inf(arb_radref(y)))
        return 0;

    if (arf_is_inf(arb_midref(x)) || arf_is_inf(arb_midref(y)) ||
        (mag_is_zero(arb_radref(x)) && mag_is_zero(arb_radref(y))))
        return arf_cmp(arb_midref(x), arb_midref(y)) < 0;

    /* xm + xr < ym - yr  <=>  xm + xr - ym + yr < 0 */
    arf_init_set_shallow(u + 0, arb_midref(x));
    arf_init_neg_shallow(u + 1, arb_midref(y));
    arf_init_set_mag_shallow(u + 2, arb_radref(x));
    arf_init_set_mag_shallow(u + 3, arb_radref(y));

    arf_init(t);
    arf_sum(t, u, 4, MAG_BITS, ARF_RND_DOWN);
    res = (arf_sgn(t) < 0);
    arf_clear(t);

    return res;
}

int arb_le(const arb_t x, const arb_t y)
{
    arf_struct u[4];
    arf_t t;
    int res;

    if (arb_is_indeterminate(x) || arb_is_indeterminate(y))
        return 0;

    if (mag_is_inf(arb_radref(x)) || mag_is_inf(arb_radref(y)))
        return (arf_is_neg_inf(arb_midref(x)) && mag_is_finite(arb_radref(x)))
            || (arf_is_pos_inf(arb_midref(y)) && mag_is_finite(arb_radref(y)));

    if (arf_is_inf(arb_midref(x)) || arf_is_inf(arb_midref(y)) ||
        (mag_is_zero(arb_radref(x)) && mag_is_zero(arb_radref(y))))
        return arf_cmp(arb_midref(x), arb_midref(y)) <= 0;

    /* xm + xr <= ym - yr  <=>  xm + xr - ym + yr <= 0 */
    arf_init_set_shallow(u + 0, arb_midref(x));
    arf_init_neg_shallow(u + 1, arb_midref(y));
    arf_init_set_mag_shallow(u + 2, arb_radref(x));
    arf_init_set_mag_shallow(u + 3, arb_radref(y));

    arf_init(t);
    arf_sum(t, u, 4, MAG_BITS, ARF_RND_DOWN);
    res = (arf_sgn(t) <= 0);
    arf_clear(t);

    return res;
}

int arb_gt(const arb_t x, const arb_t y)
{
    arf_struct u[4];
    arf_t t;
    int res;

    if (arb_is_indeterminate(x) || arb_is_indeterminate(y))
        return 0;

    if (mag_is_inf(arb_radref(x)) || mag_is_inf(arb_radref(y)))
        return 0;

    if (arf_is_inf(arb_midref(x)) || arf_is_inf(arb_midref(y)) ||
        (mag_is_zero(arb_radref(x)) && mag_is_zero(arb_radref(y))))
        return arf_cmp(arb_midref(x), arb_midref(y)) > 0;

    /* xm - xr > ym + yr  <=>  xm - xr - ym - yr > 0 */
    arf_init_set_shallow(u + 0, arb_midref(x));
    arf_init_neg_shallow(u + 1, arb_midref(y));
    arf_init_neg_mag_shallow(u + 2, arb_radref(x));
    arf_init_neg_mag_shallow(u + 3, arb_radref(y));

    arf_init(t);
    arf_sum(t, u, 4, MAG_BITS, ARF_RND_DOWN);
    res = (arf_sgn(t) > 0);
    arf_clear(t);

    return res;
}

int arb_ge(const arb_t x, const arb_t y)
{
    arf_struct u[4];
    arf_t t;
    int res;

    if (arb_is_indeterminate(x) || arb_is_indeterminate(y))
        return 0;

    if (mag_is_inf(arb_radref(x)) || mag_is_inf(arb_radref(y)))
        return (arf_is_pos_inf(arb_midref(x)) && mag_is_finite(arb_radref(x)))
            || (arf_is_neg_inf(arb_midref(y)) && mag_is_finite(arb_radref(y)));

    if (arf_is_inf(arb_midref(x)) || arf_is_inf(arb_midref(y)) ||
        (mag_is_zero(arb_radref(x)) && mag_is_zero(arb_radref(y))))
        return arf_cmp(arb_midref(x), arb_midref(y)) >= 0;

    /* xm - xr >= ym + yr  <=>  xm - xr - ym - yr >= 0 */
    arf_init_set_shallow(u + 0, arb_midref(x));
    arf_init_neg_shallow(u + 1, arb_midref(y));
    arf_init_neg_mag_shallow(u + 2, arb_radref(x));
    arf_init_neg_mag_shallow(u + 3, arb_radref(y));

    arf_init(t);
    arf_sum(t, u, 4, MAG_BITS, ARF_RND_DOWN);
    res = (arf_sgn(t) >= 0);
    arf_clear(t);

    return res;
}

