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

static void
_fmpcb_gamma(fmpcb_t y, const fmpcb_t x, long prec, int inverse)
{
    int reflect;
    long r, n, wp;
    fmpcb_t t, u, v;

    wp = prec + FLINT_BIT_COUNT(prec);

    gamma_stirling_choose_param_fmpcb(&reflect, &r, &n, x, 1, 0, wp);

    fmpcb_init(t);
    fmpcb_init(u);
    fmpcb_init(v);

    if (reflect)
    {
        /* gamma(x) = (rf(1-x, r) * pi) / (gamma(1-x+r) sin(pi x)) */
        fmpcb_sub_ui(t, x, 1, wp);
        fmpcb_neg(t, t);
        gamma_rising_fmpcb_ui_bsplit(u, t, r, wp);
        fmprb_const_pi(fmpcb_realref(v), wp);
        fmpcb_mul_fmprb(u, u, fmpcb_realref(v), wp);
        fmpcb_add_ui(t, t, r, wp);
        gamma_stirling_eval_fmpcb(v, t, n, 0, wp);
        fmpcb_exp(v, v, wp);
        fmpcb_sin_pi(t, x, wp);
        fmpcb_mul(v, v, t, wp);
    }
    else
    {
        /* gamma(x) = gamma(x+r) / rf(x,r) */
        fmpcb_add_ui(t, x, r, wp);
        gamma_stirling_eval_fmpcb(u, t, n, 0, wp);
        fmpcb_exp(u, u, prec);
        gamma_rising_fmpcb_ui_bsplit(v, x, r, wp);
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

