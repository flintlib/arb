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

#include "fmpcb.h"
#include "zeta.h"

void
fmpcb_hurwitz_zeta(fmpcb_t z, const fmpcb_t s, const fmpcb_t a, long prec)
{
    zeta_series(z, s, a, 0, 1, prec);
}

void
fmpcb_zeta(fmpcb_t z, const fmpcb_t s, long prec)
{
    fmpcb_t a;
    fmpcb_init(a);
    fmpcb_one(a);

    if (fmpr_sgn(fmprb_midref(fmpcb_realref(s))) < 0)
    {
        fmpcb_t t, u, v;
        long wp = prec + 6;

        fmpcb_init(t);
        fmpcb_init(u);
        fmpcb_init(v);

        fmpcb_sub_ui(t, s, 1, wp);

        /* 2 * (2pi)^(s-1) */
        fmprb_const_pi(fmpcb_realref(u), wp);
        fmpcb_mul_2exp_si(u, u, 1);
        fmpcb_pow(u, u, t, wp);
        fmpcb_mul_2exp_si(u, u, 1);

        /* sin(pi*s/2) */
        fmpcb_mul_2exp_si(v, s, -1);
        fmpcb_sin_pi(v, v, wp);
        fmpcb_mul(u, u, v, wp);

        /* gamma(1-s) zeta(1-s) */
        fmpcb_neg(t, t);
        fmpcb_gamma(v, t, wp);
        fmpcb_mul(u, u, v, wp);
        fmpcb_hurwitz_zeta(v, t, a, wp);
        fmpcb_mul(z, u, v, prec);

        fmpcb_clear(t);
        fmpcb_clear(u);
        fmpcb_clear(v);
    }
    else
    {
        fmpcb_hurwitz_zeta(z, s, a, prec);
    }

    fmpcb_clear(a);
}

