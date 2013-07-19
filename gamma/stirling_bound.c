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

#include <math.h>
#include "gamma.h"
#include "bernoulli.h"

void fmpr_gamma_ui_lbound(fmpr_t x, ulong n, long prec);

void fmpr_gamma_ui_ubound(fmpr_t x, ulong n, long prec);

/*
  2 |B_{2n}| G(2n+k-1) / (G(k+1) G(2n+1)) |z| (T |z|^{-1})^(2n+k)
*/
void
gamma_stirling_bound_fmpcb(fmpr_struct * err, const fmpcb_t z, long k0, long knum, long n)
{
    fmpr_t c, t, u, v;
    long i, k, prec = FMPRB_RAD_PREC;

    if (fmprb_contains_zero(fmpcb_imagref(z)) &&
        fmprb_contains_nonpositive(fmpcb_realref(z)))
    {
        for (i = 0; i < knum; i++)
            fmpr_pos_inf(err + i);
        return;
    }

    fmpr_init(c);
    fmpr_init(t);
    fmpr_init(u);
    fmpr_init(v);

    /* t = lower bound for |z| */
    fmpcb_get_abs_lbound_fmpr(t, z, prec);
    fmpcb_get_abs_ubound_fmpr(v, z, prec);

    /* c = upper bound for 1/(cos(arg(z)/2) |z|) */
    gamma_stirling_bound_phase(c, z, prec);
    fmpr_div(c, c, t, prec, FMPR_RND_UP);

    /* numerator: 2 B_{2n} gamma(2n+k-1) |z| */
    BERNOULLI_ENSURE_CACHED(2 * n);
    fmpr_set_round_fmpz(err, fmpq_numref(bernoulli_cache + 2 * n), prec, FMPR_RND_UP);
    fmpr_abs(err, err);
    fmpr_div_fmpz(err, err, fmpq_denref(bernoulli_cache + 2 * n), prec, FMPR_RND_UP);
    fmpr_mul_2exp_si(err, err, 1);
    fmpr_gamma_ui_ubound(u, 2 * n + k0 - 1, prec);
    fmpr_mul(err, err, u, prec, FMPR_RND_UP);
    fmpr_mul(err, err, v, prec, FMPR_RND_UP);

    /* denominator gamma(k+1) gamma(2n+1) */
    fmpr_gamma_ui_lbound(t, 2 * n + 1, prec);
    fmpr_gamma_ui_lbound(u, k0 + 1, prec);
    fmpr_mul(t, t, u, prec, FMPR_RND_DOWN);
    fmpr_div(err, err, t, prec, FMPR_RND_UP);

    /* multiply by c^(2n+k) */
    fmpr_pow_sloppy_ui(t, c, 2 * n + k0, prec, FMPR_RND_UP);
    fmpr_mul(err, err, t, prec, FMPR_RND_UP);

    for (i = 1; i < knum; i++)
    {
        /* recurrence factor: c * (2n+k-2) / k */
        k = k0 + i;
        fmpr_mul(err + i, err + i - 1, c, prec, FMPR_RND_UP);
        fmpr_mul_ui(err + i, err + i, 2 * n + k - 2, prec, FMPR_RND_UP);
        fmpr_div_ui(err + i, err + i, k, prec, FMPR_RND_UP);
    }

    fmpr_clear(c);
    fmpr_clear(t);
    fmpr_clear(u);
    fmpr_clear(v);
}

void
gamma_stirling_bound_fmprb(fmpr_struct * err, const fmprb_t x, long k0, long knum, long n)
{
    fmpcb_t z;
    fmpcb_init(z);
    fmpcb_set_fmprb(z, x);
    gamma_stirling_bound_fmpcb(err, z, k0, knum, n);
    fmpcb_clear(z);
}

