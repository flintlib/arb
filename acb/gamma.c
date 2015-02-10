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

#include "bernoulli.h"
#include "acb.h"

void
acb_gamma_stirling_choose_param(int * reflect, long * r, long * n,
    const acb_t z, int use_reflect, int digamma, long prec);

void acb_gamma_stirling_bound(mag_ptr err, const acb_t z, long k0, long knum, long n);

void arb_gamma_stirling_coeff(arb_t b, ulong k, int digamma, long prec);

void
acb_gamma_stirling_eval(acb_t s, const acb_t z, long nterms, int digamma, long prec)
{
    acb_t t, logz, zinv, zinv2;
    arb_t b;
    mag_t err;

    long k, term_prec;
    double z_mag, term_mag;

    acb_init(t);
    acb_init(logz);
    acb_init(zinv);
    acb_init(zinv2);
    arb_init(b);

    acb_log(logz, z, prec);
    acb_inv(zinv, z, prec);

    nterms = FLINT_MAX(nterms, 1);

    acb_zero(s);
    if (nterms > 1)
    {
        acb_mul(zinv2, zinv, zinv, prec);

        z_mag = arf_get_d(arb_midref(acb_realref(logz)), ARF_RND_UP) * 1.44269504088896;

        for (k = nterms - 1; k >= 1; k--)
        {
            term_mag = bernoulli_bound_2exp_si(2 * k);
            term_mag -= (2 * k - 1) * z_mag;
            term_prec = prec + term_mag;
            term_prec = FLINT_MIN(term_prec, prec);
            term_prec = FLINT_MAX(term_prec, 10);

            arb_gamma_stirling_coeff(b, k, digamma, term_prec);

            if (prec > 2000)
            {
                acb_set_round(t, zinv2, term_prec);
                acb_mul(s, s, t, term_prec);
            }
            else
                acb_mul(s, s, zinv2, term_prec);

            arb_add(acb_realref(s), acb_realref(s), b, term_prec);
        }

        if (digamma)
            acb_mul(s, s, zinv2, prec);
        else
            acb_mul(s, s, zinv, prec);
    }

    /* remainder bound */
    mag_init(err);
    acb_gamma_stirling_bound(err, z, digamma ? 1 : 0, 1, nterms);
    mag_add(arb_radref(acb_realref(s)), arb_radref(acb_realref(s)), err);
    mag_add(arb_radref(acb_imagref(s)), arb_radref(acb_imagref(s)), err);
    mag_clear(err);

    if (digamma)
    {
        acb_neg(s, s);
        acb_mul_2exp_si(zinv, zinv, -1);
        acb_sub(s, s, zinv, prec);
        acb_add(s, s, logz, prec);
    }
    else
    {
        /* (z-0.5)*log(z) - z + log(2*pi)/2 */
        arb_one(b);
        arb_mul_2exp_si(b, b, -1);
        arb_set(acb_imagref(t), acb_imagref(z));
        arb_sub(acb_realref(t), acb_realref(z), b, prec);
        acb_mul(t, logz, t, prec);
        acb_add(s, s, t, prec);
        acb_sub(s, s, z, prec);
        arb_const_log_sqrt2pi(b, prec);
        arb_add(acb_realref(s), acb_realref(s), b, prec);
    }

    acb_clear(t);
    acb_clear(logz);
    acb_clear(zinv);
    acb_clear(zinv2);
    arb_clear(b);
}

static void
_acb_gamma(acb_t y, const acb_t x, long prec, int inverse)
{
    int reflect;
    long r, n, wp;
    acb_t t, u, v;

    wp = prec + FLINT_BIT_COUNT(prec);

    acb_gamma_stirling_choose_param(&reflect, &r, &n, x, 1, 0, wp);

    acb_init(t);
    acb_init(u);
    acb_init(v);

    if (reflect)
    {
        /* gamma(x) = (rf(1-x, r) * pi) / (gamma(1-x+r) sin(pi x)) */
        acb_sub_ui(t, x, 1, wp);
        acb_neg(t, t);
        acb_rising_ui_rec(u, t, r, wp);
        arb_const_pi(acb_realref(v), wp);
        acb_mul_arb(u, u, acb_realref(v), wp);
        acb_add_ui(t, t, r, wp);
        acb_gamma_stirling_eval(v, t, n, 0, wp);
        acb_exp(v, v, wp);
        acb_sin_pi(t, x, wp);
        acb_mul(v, v, t, wp);
    }
    else
    {
        /* gamma(x) = gamma(x+r) / rf(x,r) */
        acb_add_ui(t, x, r, wp);
        acb_gamma_stirling_eval(u, t, n, 0, wp);
        acb_exp(u, u, prec);
        acb_rising_ui_rec(v, x, r, wp);
    }

    if (inverse)
        acb_div(y, v, u, prec);
    else
        acb_div(y, u, v, prec);

    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
}

void
acb_gamma(acb_t y, const acb_t x, long prec)
{
    if (acb_is_real(x))
    {
        arb_gamma(acb_realref(y), acb_realref(x), prec);
        arb_zero(acb_imagref(y));
        return;
    }

    _acb_gamma(y, x, prec, 0);
}

void
acb_rgamma(acb_t y, const acb_t x, long prec)
{
    if (acb_is_real(x))
    {
        arb_rgamma(acb_realref(y), acb_realref(x), prec);
        arb_zero(acb_imagref(y));
        return;
    }

    _acb_gamma(y, x, prec, 1);
}

/* corrects branch cut of sum_{k=0}^{r-1} log(z+k), given the
   logarithm of the product */
void
_acb_log_rising_correct_branch(acb_t t,
        const acb_t t_wrong, const acb_t z, ulong r, long prec)
{
    acb_t f;
    arb_t pi, u, v;
    fmpz_t pi_mult;
    long i, argprec;

    acb_init(f);

    arb_init(u);
    arb_init(pi);
    arb_init(v);

    fmpz_init(pi_mult);

    argprec = FLINT_MIN(prec, 40);

    arb_zero(u);
    for (i = 0; i < r; i++)
    {
        acb_add_ui(f, z, i, argprec);
        acb_arg(v, f, argprec);
        arb_add(u, u, v, argprec);
    }

    if (argprec == prec)
    {
        arb_set(acb_imagref(t), u);
    }
    else
    {
        arb_sub(v, u, acb_imagref(t), argprec);
        arb_const_pi(pi, argprec);
        arb_div(v, v, pi, argprec);

        if (arb_get_unique_fmpz(pi_mult, v))
        {
            arb_const_pi(v, prec);
            arb_mul_fmpz(v, v, pi_mult, prec);
            arb_add(acb_imagref(t), acb_imagref(t), v, prec);
        }
        else
        {
            arb_zero(u);
            for (i = 0; i < r; i++)
            {
                acb_add_ui(f, z, i, prec);
                acb_arg(v, f, prec);
                arb_add(u, u, v, prec);
            }
            arb_set(acb_imagref(t), u);
        }
    }

    acb_clear(f);

    arb_clear(u);
    arb_clear(v);
    arb_clear(pi);

    fmpz_clear(pi_mult);
}

void
acb_lgamma(acb_t y, const acb_t x, long prec)
{
    int reflect;
    long r, n, wp;
    acb_t t, u;

    if (acb_is_real(x) && arb_is_positive(acb_realref(x)))
    {
        arb_lgamma(acb_realref(y), acb_realref(x), prec);
        arb_zero(acb_imagref(y));
        return;
    }

    wp = prec + FLINT_BIT_COUNT(prec);

    acb_gamma_stirling_choose_param(&reflect, &r, &n, x, 0, 0, wp);

    /* log(gamma(x)) = log(gamma(x+r)) - log(rf(x,r)) */
    acb_init(t);
    acb_init(u);

    acb_add_ui(t, x, r, wp);
    acb_gamma_stirling_eval(u, t, n, 0, wp);

    acb_rising_ui_rec(t, x, r, prec);
    acb_log(t, t, prec);

    _acb_log_rising_correct_branch(t, t, x, r, wp);

    acb_sub(y, u, t, prec);

    acb_clear(t);
    acb_clear(u);
}

