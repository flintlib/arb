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

/* |B_n| / n! < 4 * (2 pi)^(-n) = 2^(2-log_2(2pi)*n) */
void
fmpr_bernoulli_scaled_ui_ubound(fmpr_t b, ulong n, long prec)
{
    fmpr_set_si_2exp_si(b, 1, 2 - 265*n/100);
}

/* Absolute value of rising factorial (could speed up once complex gamma is available). */
void
fmpcb_rfac_abs_ubound(fmpr_t bound, const fmpcb_t s, ulong n, long prec)
{
    fmpr_t term, t;
    ulong k;

    /* M(k) = (a+k)^2 + b^2
       M(0) = a^2 + b^2
       M(k+1) = M(k) + 2*a + (2*k+1)
    */
    fmpr_init(t);
    fmpr_init(term);

    fmpr_one(bound);

    /* M(0) = a^2 + b^2 */
    fmprb_get_abs_ubound_fmpr(t, fmpcb_realref(s), prec);
    fmpr_mul(term, t, t, prec, FMPR_RND_UP);
    fmprb_get_abs_ubound_fmpr(t, fmpcb_imagref(s), prec);
    fmpr_mul(t, t, t, prec, FMPR_RND_UP);
    fmpr_add(term, term, t, prec, FMPR_RND_UP);

    /* we add t = 2*a to each term. note that this can be signed;
       we always want the most positive value */
    fmpr_add(t, fmprb_midref(fmpcb_realref(s)),
        fmprb_radref(fmpcb_realref(s)), prec, FMPR_RND_CEIL);
    fmpr_mul_2exp_si(t, t, 1);

    for (k = 0; k < n; k++)
    {
        fmpr_mul(bound, bound, term, prec, FMPR_RND_UP);
        fmpr_add_ui(term, term, 2 * k + 1, prec, FMPR_RND_UP);
        fmpr_add(term, term, t, prec, FMPR_RND_UP);
    }

    fmpr_sqrt(bound, bound, prec, FMPR_RND_UP);

    fmpr_clear(t);
    fmpr_clear(term);
}

void
_fmpr_get_floor_fmpz(fmpz_t f, const fmpr_t x)
{
    if (COEFF_IS_MPZ(*fmpr_expref(x)))
    {
        abort();
    }

    if (fmpz_sgn(fmpr_expref(x)) >= 0)
    {
        fmpz_mul_2exp(f, fmpr_manref(x), *fmpr_expref(x));
    }
    else
    {
        fmpz_neg(f, fmpr_expref(x));
        fmpz_fdiv_q_2exp(f, fmpr_manref(x), *f);
    }
}

/* the bound is rf(s,2M) * (B_{2M} / (2M)!) / ((a+2M-1) * N^(a+2M-1))
   we assume N > 0, Q > 0 */
void
fmpcb_zeta_em_bound(fmpr_t bound, const fmpcb_t s, ulong N, ulong M, long prec)
{
    fmpr_t delta, t, u;

    fmpr_init(delta);

    /* delta = a + 2*q - 1, rounded down (we require it to be positive) */
    fmpr_sub(delta, fmprb_midref(fmpcb_realref(s)),
        fmprb_radref(fmpcb_realref(s)), prec, FMPR_RND_FLOOR);
    fmpr_add_ui(delta, delta, 2*M - 1, prec, FMPR_RND_FLOOR);

    if (fmpr_sgn(delta) <= 0 || N < 1 || M < 1)
    {
        fmpr_clear(delta);
        fmpr_pos_inf(bound);
        return;
    }

    fmpr_init(t);
    fmpr_init(u);

    /* compute numerator (rounding up) */
    fmpcb_rfac_abs_ubound(t, s, 2 * M, prec);
    fmpr_bernoulli_scaled_ui_ubound(u, 2*M, prec);
    fmpr_mul(t, t, u, prec, FMPR_RND_UP);

    /* compute denominator (rounding down) */
    fmpr_set_ui(u, N);
    {
        fmpz_t f;
        fmpz_init(f);
        _fmpr_get_floor_fmpz(f, delta);
        fmpr_pow_sloppy_fmpz(u, u, f, prec, FMPR_RND_DOWN);
        fmpz_clear(f);
    }

    fmpr_mul(u, u, delta, prec, FMPR_RND_DOWN);
    fmpr_div(bound, t, u, prec, FMPR_RND_UP);

    fmpr_clear(delta);
    fmpr_clear(t);
    fmpr_clear(u);
}
