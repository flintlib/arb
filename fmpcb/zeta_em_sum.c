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
#include "bernoulli.h"

void
fmpcb_zeta_em_sum(fmpcb_t z, const fmpcb_t s, ulong N, ulong M, long prec)
{
    fmpcb_t t, u, v, negs, term, sum;
    fmprb_t x;
    fmpz_t c;
    ulong r, n;

    bernoulli_cache_compute(2 * M);

    fmpcb_init(t);
    fmpcb_init(u);
    fmpcb_init(v);
    fmpcb_init(negs);
    fmpcb_init(term);
    fmpcb_init(sum);

    fmprb_init(x);
    fmpz_init(c);

    fmpcb_neg(negs, s);
    fmpcb_zero(sum);

    /* sum 1/n^s */
    for (n = 1; n < N; n++)
    {
        fmpcb_set_ui(t, n);
        fmpcb_pow(t, t, negs, prec);
        fmpcb_add(sum, sum, t, prec);
    }

    /* t = 1/N^s */
    fmpcb_set_ui(t, N);
    fmpcb_pow(t, t, negs, prec);

    /* N / (s-1) * (1 / N^s) */
    fmpcb_set_ui(u, N);
    fmpcb_sub_ui(v, s, 1, prec);
    fmpcb_div(u, u, v, prec);
    fmpcb_mul(u, u, t, prec);
    fmpcb_add(sum, sum, u, prec);

    /* 1/2 * (1 / N^s) */
    fmpcb_mul_2exp_si(u, t, -1);
    fmpcb_add(sum, sum, u, prec);

    /* term = 1/2 * (1 / N^s) * s / N */
    fmpcb_mul(u, u, s, prec);
    fmpcb_div_ui(term, u, N, prec);

    for (r = 1; r < M; r++)
    {
        /* sum += bernoulli number * term */
        fmprb_set_round_fmpz(x, fmpq_numref(bernoulli_cache + 2 * r), prec);
        fmprb_div_fmpz(x, x, fmpq_denref(bernoulli_cache + 2 * r), prec);

        fmprb_mul(fmpcb_realref(u), fmpcb_realref(term), x, prec);
        fmprb_mul(fmpcb_imagref(u), fmpcb_imagref(term), x, prec);

        fmpcb_add(sum, sum, u, prec);

        /* multiply term by (s+2r-1)(s+2r) / (N^2 * (2*r+1)*(2*r+2)) */
        fmpcb_set(u, s);
        fmprb_add_ui(fmpcb_realref(u), fmpcb_realref(s), 2*r-1, prec);
        fmpcb_mul(term, term, u, prec);
        fmprb_add_ui(fmpcb_realref(u), fmpcb_realref(u), 1, prec);

        fmpcb_mul(term, term, u, prec);

        fmpz_set_ui(c, N);
        fmpz_mul_ui(c, c, N);
        fmpz_mul_ui(c, c, 2*r+1);
        fmpz_mul_ui(c, c, 2*r+2);
        fmprb_div_fmpz(fmpcb_realref(term), fmpcb_realref(term), c, prec);
        fmprb_div_fmpz(fmpcb_imagref(term), fmpcb_imagref(term), c, prec);

    }

    fmpcb_set(z, sum);

    fmpcb_clear(t);
    fmpcb_clear(u);
    fmpcb_clear(v);
    fmpcb_clear(negs);
    fmpcb_clear(term);
    fmpcb_clear(sum);

    fmprb_clear(x);
    fmpz_clear(c);
}
