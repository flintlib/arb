/*
    Copyright (C) 2022 Erik Postma

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

void _arb_poly_to_taylor_model(arb_ptr g, mag_t radius, arb_srcptr f, slong len, const arb_t x, slong glen,
    slong prec)
{
    /*
     * Consider a value x0 in x. Write x0 = xc + xr, with xc = arb_midref(x). Express f as an exact
     * polynomial g in xr = x0 - xc, with a bound, radius, to enclose all values that f can assume.
     *
     * Define the stages of constructing the Horner form of f, expressed in an abstract variable y,
     * as:
     *
     *   HP[len](y) = 0,
     *   HP[d](y) = f[d] + y * HP[d+1](y),    for 0 <= d < len.
     *
     * In addition to making sure that each g's radius is 0, we maintain this invariant in the loop
     * below:
     *
     *   abs(HP[i](x0) - sum(g[d] * xr^d, d = 0 .. glen-1)) <= radius.
     *
     * In order to establish this from the invariant for i+1, we use the definition of HP and write
     *
     *   abs((f[i] + x0 * HP[i+1](x0)) - (f[i] + x0 * sum(g[d] * xr^d, d = 0 .. glen-1))) <=
     *   abs(f[i]) * radius.
     *
     * We first assume that f[i] is exact. The first two terms inside the absolute value are
     * HP[i](x0); minus the rest can be expanded to
     *
     *   f[i] + (xc + xr) * sum(g[d] * xr^d, d = 0 .. glen-1) =
     *   f[i] + sum(xc * g[d] * xr^d, d = 0 .. glen-1)
     *        + sum(g[d] * xr^(d+1), d = 0 .. glen-1) =
     *   f[i] + xc * g[0]
     *        + sum((g[d-1] + xc * g[d]) * xr^d, d = 1 .. glen-1)
     *        + g[glen-1] * xr^glen.
     *
     * So we need to:
     *
     * - multiply radius by abs(x);
     * - set g[0] to f[i] + xc * g[0]
     * - set g[d] to g[d-1] + xc * g[d] for d = 1 .. glen-1
     * - account for g[glen-1] * xr^glen by bounding it appropriately:
     *   * if glen is odd, add abs(g[glen-1]) * arb_radref(x)^glen to radius;
     *   * if glen is even, add abs(g[glen-1]) * (arb_radref(x)^glen)/2 to radius and add g[glen-1]
     *     * (arb_radref(x)^glen)/2 to g[0].
     *
     * Finally, if f[i] has a nonzero radius, then we can just add it to radius.
     */
    arf_t g_glen_m_1;           /* Stores a copy of g[glen - 1]. */
    mag_t abs_x;                /* Precompute abs(x). */
    mag_t u;                    /* Scratch value. */
    mag_t rad_x_glen_or_half;   /* Precompute arb_radref(x)^glen, divided by 2 if glen is even. */
    arf_t rad_x_glen_half_arf;  /* If glen is even, same, but as an arf_t. This gets initialized
                                 * *only* if glen is even. */
    slong i, d;
    
    arf_init(g_glen_m_1);

    mag_init(abs_x);
    arb_get_mag(abs_x, x);

    mag_init(u);

    mag_init(rad_x_glen_or_half);
    mag_pow_ui(rad_x_glen_or_half, arb_radref(x), glen);
    if(glen % 2 == 0)
    {
        mag_div_ui(rad_x_glen_or_half, rad_x_glen_or_half, 2);
        arf_init(rad_x_glen_half_arf);
        arf_set_mag(rad_x_glen_half_arf, rad_x_glen_or_half);
    }
    
    _arb_vec_zero(g, glen);
    mag_zero(radius);

    for(i = len - 1; i >= 0; --i)
    {
        /* Set radius to radius * abs(x) + abs(g[glen-1]) * arb_radref(x)^glen (maybe divided by 2)
         * + arb_radref(f[i]). */
        mag_mul(radius, radius, abs_x);
        arf_get_mag(u, g_glen_m_1);
        mag_addmul(radius, u, rad_x_glen_or_half);
        mag_add(radius, radius, arb_radref(f + i));

        /* We need to overwrite g[glen-1] before we can use it to update g[0]. */
        arf_set(g_glen_m_1, arb_midref(g + (glen-1)));

        for(d = glen - 1; d >= 1; --d)
        {
            /* Set g[d] to g[d-1] + xc * g[d]. */
            arf_fma(arb_midref(g + d), arb_midref(x), arb_midref(g + d), arb_midref(g + (d-1)),
                prec, ARF_RND_NEAR);
        }

        /* Set g[0] to f[i] + xc * g[0] (maybe plus g[glen-1] * arb_radref(x)^glen / 2). */
        arf_fma(arb_midref(g + 0), arb_midref(x), arb_midref(g + 0), arb_midref(f + i),
            prec, ARF_RND_NEAR);

        if(glen % 2 == 0)
        {
            arf_addmul(arb_midref(g + 0), g_glen_m_1, rad_x_glen_half_arf, prec, ARF_RND_NEAR);
        }

        /*
        for(d = 0; d < glen; ++d)
        {
            flint_printf("g[%d] = ", d);   arb_printd(g + d, 15);  flint_printf("\n");
        }
        flint_printf("radius = ");         mag_printd(radius, 15); flint_printf("\n\n");
        */
    }

    if(glen % 2 == 0)
    {
        arf_clear(rad_x_glen_half_arf);
    }
    mag_clear(rad_x_glen_or_half);
    mag_clear(u);
    mag_clear(abs_x);
    arf_clear(g_glen_m_1);
}
