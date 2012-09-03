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
fmprb_bin_ui_bsplit(fmprb_t P, fmprb_t Q,
    const fmprb_t n, ulong a, ulong b, ulong k, long prec)
{
    if (b - a == 1)
    {
        fmprb_sub_ui(P, n, k - b, prec);
        fmprb_set_ui(Q, b);
    }
    else
    {
        fmprb_t R, S;
        long m = a + (b - a) / 2;

        fmprb_init(R);
        fmprb_init(S);
        fmprb_bin_ui_bsplit(P, Q, n, a, m, k, prec);
        fmprb_bin_ui_bsplit(R, S, n, m, b, k, prec);
        fmprb_mul(P, P, R, prec);
        fmprb_mul(Q, Q, S, prec);

        fmprb_clear(R);
        fmprb_clear(S);
    }
}

void
fmprb_bin_ui(fmprb_t x, const fmprb_t n, ulong k, long prec)
{
    if (k == 0)
    {
        fmprb_set_ui(x, 1UL);
    }
    else
    {
        fmprb_t P, Q;

        fmprb_init(P);
        fmprb_init(Q);

        fmprb_bin_ui_bsplit(P, Q, n, 0, k, k, prec + FLINT_BIT_COUNT(k));

        fmprb_div(x, P, Q, prec);

        fmprb_clear(P);
        fmprb_clear(Q);
    }
}

void
fmprb_bin_uiui(fmprb_t x, ulong n, ulong k, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_set_ui(t, n);
    fmprb_bin_ui(x, t, k, prec);
    fmprb_clear(t);
}

