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

#include "fmpcb.h"
#include "gamma.h"

/* tuning factor */
#define GAMMA_STIRLING_BETA 0.21

static void
choose(int * reflect, long * r, long * n, const fmpcb_t z, int use_refl, long prec)
{
    double x, y;

    x = fmpr_get_d(fmprb_midref(fmpcb_realref(z)), FMPR_RND_NEAR);
    y = fmpr_get_d(fmprb_midref(fmpcb_imagref(z)), FMPR_RND_NEAR);

    gamma_stirling_choose_param(reflect, r, n, x, y,
        GAMMA_STIRLING_BETA, use_refl, prec);
}

static void
_fmpcb_gamma(fmpcb_t y, const fmpcb_t x, long prec, int inverse)
{
    int reflect;
    long r, n, wp;
    fmpcb_t t, u, v;

    wp = prec + FLINT_BIT_COUNT(prec);

    choose(&reflect, &r, &n, x, 1, wp);

    fmpcb_init(t);
    fmpcb_init(u);
    fmpcb_init(v);

    if (reflect)
    {
        /* gamma(x) = (rf(1-x, r) * pi) / (gamma(1-x+r) sin(pi x)) */
        fmpcb_sub_ui(t, x, 1, wp);
        fmpcb_neg(t, t);
        fmpcb_rfac_ui_bsplit(u, t, r, wp);
        fmprb_const_pi(fmpcb_realref(v), wp);
        fmpcb_mul_fmprb(u, u, fmpcb_realref(v), wp);
        fmpcb_add_ui(t, t, r, wp);
        gamma_stirling_eval_series_fmpcb(v, t, n, wp);
        fmpcb_exp(v, v, wp);
        fmpcb_sin_pi(t, x, wp);
        fmpcb_mul(v, v, t, wp);
    }
    else
    {
        /* gamma(x) = gamma(x+r) / rf(x,r) */
        fmpcb_add_ui(t, x, r, wp);
        gamma_stirling_eval_series_fmpcb(u, t, n, wp);
        fmpcb_exp(u, u, prec);
        fmpcb_rfac_ui_bsplit(v, x, r, wp);
    }

    if (inverse)
        fmpcb_div(y, v, u, prec);
    else
        fmpcb_div(y, u, v, prec);

    fmpcb_clear(t);
    fmpcb_clear(u);
    fmpcb_clear(v);
}

void
fmpcb_gamma(fmpcb_t y, const fmpcb_t x, long prec)
{
    _fmpcb_gamma(y, x, prec, 0);
}

void
fmpcb_rgamma(fmpcb_t y, const fmpcb_t x, long prec)
{
    _fmpcb_gamma(y, x, prec, 1);
}

void
fmpcb_log_rfac_ui(fmpcb_t s, const fmpcb_t z, ulong r, long wp)
{
    ulong i;
    fmpcb_t t, u;

    fmpcb_init(t);
    fmpcb_init(u);

    fmpcb_zero(t);

    for (i = 0; i < r; i++)
    {
        fmpcb_add_ui(u, z, i, wp);
        fmpcb_log(u, u, wp);
        fmpcb_add(t, t, u, wp);
    }

    fmpcb_set(s, t);

    fmpcb_clear(t);
    fmpcb_clear(u);
}

void
fmpcb_lgamma(fmpcb_t y, const fmpcb_t x, long prec)
{
    int reflect;
    long r, n, wp;
    fmpcb_t t, u;

    wp = prec + FLINT_BIT_COUNT(prec);

    choose(&reflect, &r, &n, x, 0, wp);

    /* log(gamma(x)) = log(gamma(x+r)) - log(rf(x,r)) */
    fmpcb_init(t);
    fmpcb_init(u);

    fmpcb_add_ui(t, x, r, wp);
    gamma_stirling_eval_series_fmpcb(u, t, n, wp);

    fmpcb_log_rfac_ui(t, x, r, wp);
    fmpcb_sub(y, u, t, prec);

    fmpcb_clear(t);
    fmpcb_clear(u);
}

