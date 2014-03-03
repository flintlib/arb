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

#include "partitions.h"

static void
_fmpr_sinh(fmpr_t y, const fmpr_t x, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_set_fmpr(t, x);
    fmprb_sinh(t, t, prec);
    fmpr_add(y, fmprb_midref(t), fmprb_radref(t), prec, FMPR_RND_UP);
    fmprb_clear(t);
}

/* Equation (1.8) in the paper */
void
partitions_rademacher_bound(fmpr_t b, const fmpz_t n, ulong N)
{
    fmpr_t A, B, C, t, u;
    fmpz_t n1;

    fmpr_init(A);
    fmpr_init(B);
    fmpr_init(C);
    fmpr_init(t);
    fmpr_init(u);
    fmpz_init(n1);

    /* bound for 44*pi^2/(225*sqrt(3)) */
    fmpr_set_si_2exp_si(A, 18695160, -24);

    /* bound for pi*sqrt(2)/75 */
    fmpr_set_si_2exp_si(B, 993857, -24);

    /* bound for pi*sqrt(2/3) */
    fmpr_set_si_2exp_si(C, 43035232, -24);

    /* first term: A / sqrt(N) */
    fmpr_sqrt_ui(t, N, FMPRB_RAD_PREC, FMPR_RND_DOWN);
    fmpr_div(b, A, t, FMPRB_RAD_PREC, FMPR_RND_UP);

    /* B * sqrt(N/(n-1)) */
    fmpr_set_ui(t, N);
    fmpz_sub_ui(n1, n, 1);
    fmpr_div_fmpz(t, t, n1, FMPRB_RAD_PREC, FMPR_RND_UP);
    fmpr_sqrt(t, t, FMPRB_RAD_PREC, FMPR_RND_UP);
    fmpr_mul(t, B, t, FMPRB_RAD_PREC, FMPR_RND_UP);

    /* sinh(C*sqrt(n)/N) */
    fmpr_sqrt_fmpz(u, n, FMPRB_RAD_PREC, FMPR_RND_UP);
    fmpr_div_ui(u, u, N, FMPRB_RAD_PREC, FMPR_RND_UP);
    fmpr_mul(u, C, u, FMPRB_RAD_PREC, FMPR_RND_UP);

    _fmpr_sinh(u, u, FMPRB_RAD_PREC);

    /* second term: B * ... * sinh... */
    fmpr_mul(t, t, u, FMPRB_RAD_PREC, FMPR_RND_UP);
    fmpr_add(b, b, t, FMPRB_RAD_PREC, FMPR_RND_UP);

    fmpr_clear(A);
    fmpr_clear(B);
    fmpr_clear(C);
    fmpr_clear(t);
    fmpr_clear(u);
    fmpz_clear(n1);
}

