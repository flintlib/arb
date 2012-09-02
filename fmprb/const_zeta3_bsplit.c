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

#include "fmprb.h"

static void
zeta3_bsplit(fmprb_t P, fmprb_t Q, fmprb_t T, long a, long b, long wp, int cont)
{
    if (b - a == 1)
    {
        fmprb_ui_pow_ui(P, a ? a : 1, 5, wp);

        fmprb_ui_pow_ui(Q, 2*a + 1, 5, wp);
        fmprb_mul_ui(Q, Q, a ? 32 : 64, wp);

        fmprb_set_ui(T, 205*a*a + 250*a + 77); /* xxx: fix overflow */
        if (a % 2)
            fmprb_neg(T, T);
        fmprb_mul(T, T, P, wp);
    }
    else
    {
        long m;
        fmprb_t P2, Q2, T2;

        m = (a + b) / 2;

        fmprb_init(P2);
        fmprb_init(Q2);
        fmprb_init(T2);

        zeta3_bsplit(P, Q, T, a, m, wp, 1);
        zeta3_bsplit(P2, Q2, T2, m, b, wp, 1);

        fmprb_mul(T, T, Q2, wp);
        fmprb_addmul(T, P, T2, wp);
        fmprb_mul(Q, Q, Q2, wp);
        if (cont)
            fmprb_mul(P, P, P2, wp);

        fmprb_clear(P2);
        fmprb_clear(Q2);
        fmprb_clear(T2);
    }
}

void
fmprb_const_zeta3_bsplit(fmprb_t x, long prec)
{
    long wp, N;
    fmprb_t P, Q, T;

    wp = prec + 10 + FLINT_BIT_COUNT(prec);
    N = prec / 10 + 2;

    fmprb_init(P);
    fmprb_init(Q);
    fmprb_init(T);

    zeta3_bsplit(P, Q, T, 0, N, wp, 0);
    fmprb_div(x, T, Q, prec);

    fmprb_add_error_2exp_si(x, -10 * (N - 1));

    fmprb_clear(P);
    fmprb_clear(Q);
    fmprb_clear(T);
}
