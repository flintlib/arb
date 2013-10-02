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
    long i, n, depth;

    fmprb_t pi, theta, v;
    fmpcb_t t, u;
    fmprb_t b;

    fmprb_init(pi);
    fmprb_init(theta);
    fmpcb_init(t);
    fmpcb_init(u);
    fmprb_init(v);
    fmprb_init(b);

    fmprb_const_pi(pi, prec);

    fmprb_zero_pm_inf(b);

    for (depth = 0, n = 16; depth < maxdepth; n *= 2, depth++)
    {
        fmprb_zero(b);

        for (i = 0; i < n; i++)
        {
            /* theta = 2 pi ([i,i+1]/n) */
            fmpr_set_si(fmprb_midref(theta), 2 * i + 1);
            fmpr_set_si(fmprb_radref(theta), 2);
            fmprb_div_ui(theta, theta, n, prec);
            fmprb_mul(theta, theta, pi, prec);

            fmprb_sin_cos(fmpcb_imagref(t), fmpcb_realref(t), theta, prec);
            fmpcb_mul_fmprb(t, t, radius, prec);
            fmpcb_add(t, t, x, prec);

            func(u, t, param, 1, prec);
            fmpcb_abs(v, u, prec);

            fmprb_add(b, b, v, prec);
        }

        fmprb_div_ui(b, b, n, prec);
        fmprb_div(b, b, radius, prec);

        if (fmprb_is_exact(b) || fmpr_cmp(fmprb_radref(b), fmprb_midref(b)) < 0)
            break;
    }

    fmprb_set(bound, b);

    fmprb_clear(pi);
    fmprb_clear(theta);
    fmpcb_clear(t);
    fmpcb_clear(u);
    fmprb_clear(v);
    fmprb_clear(b);
}

