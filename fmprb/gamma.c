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
#include "arith.h"
#include "fmprb.h"
#include "bernoulli.h"

void
fmprb_stirling_series_coeff(fmprb_t b, ulong k, long prec)
{
    fmpz_t d;
    fmpz_init(d);
    fmprb_set_round_fmpz(b, fmpq_numref(bernoulli_cache + 2 * k), prec);
    fmpz_mul2_uiui(d, fmpq_denref(bernoulli_cache + 2 * k), 2 * k, 2 * k - 1);
    fmprb_div_fmpz(b, b, d, prec);
    fmpz_clear(d);
}

/*
Heuristic: we use Stirling's series if abs(x) > beta * prec.

For convergence, beta must be greater than log(2)/(2*pi) ~= 0.11.
A larger beta gives faster convergence at the expense of extra
argument reduction.
*/

#define GAMMA_STIRLING_BETA 0.21

long stirling_choose_r(const fmprb_t x, long wp)
{
    double t = fmpr_get_d(fmprb_midref(x), FMPR_RND_FLOOR);
    double want = FLINT_MAX(5, GAMMA_STIRLING_BETA * wp);

    return (long) FLINT_MAX(0, want - t + 1);
}

/* TODO: speed up by caching Bernoulli number sizes */
long
stirling_choose_nterms(const fmprb_t x, long r, double bits)
{
    long i;
    double t, logt;
    double mag;

    t = fmpr_get_d(fmprb_midref(x), FMPR_RND_FLOOR) + r;
    logt = log(t);

    for (i = 1; ; i++)
    {
        mag = bernoulli_bound_2exp_si(2 * i) - (logt * 1.4426950408889634) * (2 * i - 1);
        if (mag < -bits)
            return i;
    }
}

void
fmprb_gamma_log_stirling(fmprb_t s, const fmprb_t z, long nterms, long prec)
{
    fmprb_t t, u, b, w, v;
    long k, term_prec;
    double z_mag, term_mag;

    fmprb_init(t);
    fmprb_init(u);
    fmprb_init(b);
    fmprb_init(w);
    fmprb_init(v);

    fmprb_log(w, z, prec);

    bernoulli_cache_compute(2 * (nterms + 1));

    nterms = FLINT_MAX(nterms, 1);

    fmprb_zero(s);
    if (nterms > 1)
    {
        fmprb_ui_div(t, 1UL, z, prec);
        fmprb_mul(u, t, t, prec);
        z_mag = fmpr_get_d(fmprb_midref(w), FMPR_RND_UP) * 1.44269504088896;

        for (k = nterms - 1; k >= 1; k--)
        {
            term_mag = bernoulli_bound_2exp_si(2 * k);
            term_mag -= (2 * k - 1) * z_mag;
            term_prec = prec + term_mag;
            term_prec = FLINT_MIN(term_prec, prec);
            term_prec = FLINT_MAX(term_prec, 10);

            fmprb_stirling_series_coeff(b, k, term_prec);

            if (prec > 2000)
            {
                fmprb_set_round(v, u, term_prec);
                fmprb_mul(s, s, v, term_prec);
            }
            else
            {
                fmprb_mul(s, s, u, term_prec);
            }

            fmprb_add(s, s, b, term_prec);
        }

        fmprb_mul(s, s, t, prec);
    }

    /* finally, we add the remainder error bound */
    /* B_2n / (2n (2n-1)) */
    fmprb_stirling_series_coeff(b, nterms, FMPRB_RAD_PREC);
    /* 1/z^(2n-1) */
    fmprb_ui_div(t, 1UL, z, FMPRB_RAD_PREC);
    fmprb_pow_ui(t, t, 2 * nterms - 1, FMPRB_RAD_PREC);
    fmprb_mul(t, t, b, FMPRB_RAD_PREC);

    fmprb_add_error(s, t);

    /* (z-0.5)*log(z) - z + log(2*pi)/2 */
    fmprb_set_ui(t, 1);
    fmprb_mul_2exp_si(t, t, -1);
    fmprb_sub(t, z, t, prec);
    fmprb_mul(t, w, t, prec);

    fmprb_add(s, s, t, prec);
    fmprb_sub(s, s, z, prec);

    fmprb_const_log_sqrt2pi(t, prec);
    fmprb_add(s, s, t, prec);

    fmprb_clear(t);
    fmprb_clear(u);
    fmprb_clear(b);
    fmprb_clear(w);
    fmprb_clear(v);
}

void
fmprb_lgamma(fmprb_t y, const fmprb_t x, long prec)
{
    long r, n, wp;
    fmprb_t t, u;

    wp = prec + FLINT_BIT_COUNT(prec);

    r = stirling_choose_r(x, wp);
    n = stirling_choose_nterms(x, r, wp);

    /* log(gamma(x)) = log(gamma(x+r)) - log(rf(x,r)) */
    fmprb_init(t);
    fmprb_init(u);

    fmprb_add_ui(t, x, r, wp);
    fmprb_gamma_log_stirling(u, t, n, wp);

    fmprb_rfac_ui_bsplit(t, x, r, wp);
    fmprb_log(t, t, wp);
    fmprb_sub(y, u, t, prec);

    fmprb_clear(t);
    fmprb_clear(u);
}

void
fmprb_gamma(fmprb_t y, const fmprb_t x, long prec)
{
    long r, n, wp;
    fmprb_t t, u;

    wp = prec + FLINT_BIT_COUNT(prec);

    r = stirling_choose_r(x, wp);
    n = stirling_choose_nterms(x, r, wp);

    /* gamma(x) = gamma(x+r) / rf(x,r) */
    fmprb_init(t);
    fmprb_init(u);

    fmprb_add_ui(t, x, r, wp);
    fmprb_gamma_log_stirling(u, t, n, wp);
    fmprb_exp(u, u, prec);

    fmprb_rfac_ui_bsplit(t, x, r, wp);
    fmprb_div(y, u, t, prec);

    fmprb_clear(t);
    fmprb_clear(u);
}

void
fmprb_rgamma(fmprb_t y, const fmprb_t x, long prec)
{
    long r, n, wp;
    fmprb_t t, u;

    wp = prec + FLINT_BIT_COUNT(prec);

    r = stirling_choose_r(x, wp);
    n = stirling_choose_nterms(x, r, wp);

    /* 1/gamma(x) = rf(x,r) / gamma(x+r) */
    fmprb_init(t);
    fmprb_init(u);

    fmprb_add_ui(t, x, r, wp);
    fmprb_gamma_log_stirling(u, t, n, wp);
    fmprb_neg(u, u);
    fmprb_exp(u, u, prec);

    fmprb_rfac_ui_bsplit(t, x, r, wp);
    fmprb_mul(y, u, t, prec);

    fmprb_clear(t);
    fmprb_clear(u);
}
