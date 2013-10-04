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

#include "fmpcb_calc.h"

void
fmpcb_calc_cauchy_bound(fmprb_t bound, fmpcb_calc_func_t func, void * param,
    const fmpcb_t x, const fmprb_t radius, long maxdepth, long prec)
{
    long i, n, depth, wp;

    fmprb_t pi, theta, v, s1, c1, s2, c2, st, ct;
    fmpcb_t t, u;
    fmprb_t b;

    fmprb_init(pi);
    fmprb_init(theta);
    fmprb_init(v);

    fmprb_init(s1);
    fmprb_init(c1);
    fmprb_init(s2);
    fmprb_init(c2);
    fmprb_init(st);
    fmprb_init(ct);

    fmpcb_init(t);
    fmpcb_init(u);
    fmprb_init(b);

    wp = prec + 20;

    fmprb_const_pi(pi, wp);
    fmprb_zero_pm_inf(b);

    for (depth = 0, n = 16; depth < maxdepth; n *= 2, depth++)
    {
        fmprb_zero(b);

        /* theta = 2 pi / n */
        fmprb_div_ui(theta, pi, n, wp);
        fmprb_mul_2exp_si(theta, theta, 1);

        /* sine and cosine of i*theta and (i+1)*theta */
        fmprb_zero(s1);
        fmprb_one(c1);
        fmprb_sin_cos(st, ct, theta, wp);
        fmprb_set(s2, st);
        fmprb_set(c2, ct);

        for (i = 0; i < n; i++)
        {
            /* sine and cosine of 2 pi ([i,i+1]/n) */

            /* since we use power of two subdivision points, the
               sine and cosine are monotone on each subinterval */
            fmprb_union(fmpcb_realref(t), c1, c2, wp);
            fmprb_union(fmpcb_imagref(t), s1, s2, wp);
            fmpcb_mul_fmprb(t, t, radius, wp);
            fmpcb_add(t, t, x, prec);

            /* next angle */
            fmprb_mul(v, c2, ct, wp);
            fmprb_mul(c1, s2, st, wp);
            fmprb_sub(c1, v, c1, wp);
            fmprb_mul(v, c2, st, wp);
            fmprb_mul(s1, s2, ct, wp);
            fmprb_add(s1, v, s1, wp);
            fmprb_swap(c1, c2);
            fmprb_swap(s1, s2);

            func(u, t, param, 1, prec);
            fmpcb_abs(v, u, prec);
            fmprb_add(b, b, v, prec);
        }

        fmprb_div_ui(b, b, n, prec);

        if (fmprb_is_exact(b) || fmpr_cmp(fmprb_radref(b), fmprb_midref(b)) < 0)
            break;
    }

    fmprb_set(bound, b);

    fmprb_clear(pi);
    fmprb_clear(theta);
    fmprb_clear(v);

    fmpcb_clear(t);
    fmpcb_clear(u);
    fmprb_clear(b);

    fmprb_clear(s1);
    fmprb_clear(c1);
    fmprb_clear(s2);
    fmprb_clear(c2);
    fmprb_clear(st);
    fmprb_clear(ct);
}

