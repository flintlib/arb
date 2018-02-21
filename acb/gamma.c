/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "bernoulli.h"
#include "acb.h"

void
acb_gamma_stirling_choose_param(int * reflect, slong * r, slong * n,
    const acb_t z, int use_reflect, int digamma, slong prec);

void acb_gamma_stirling_bound(mag_ptr err, const acb_t z, slong k0, slong knum, slong n);

void arb_gamma_stirling_coeff(arb_t b, ulong k, int digamma, slong prec);

void
acb_gamma_stirling_eval(acb_t s, const acb_t z, slong nterms, int digamma, slong prec)
{
    acb_t t, logz, zinv, zinv2;
    arb_t b;
    mag_t err;

    slong k, term_prec;
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
_acb_gamma(acb_t y, const acb_t x, slong prec, int inverse)
{
    int reflect;
    slong r, n, wp;
    acb_t t, u, v;
    double acc;

    wp = prec + FLINT_BIT_COUNT(prec);

    /* todo: for large x (if exact or accurate enough), increase precision */
    acc = acb_rel_accuracy_bits(x);
    acc = FLINT_MAX(acc, 0);
    wp = FLINT_MIN(prec, acc + 20);
    wp = FLINT_MAX(wp, 2);
    wp = wp + FLINT_BIT_COUNT(wp);

    if (acc < 3)  /* try to avoid divisions blowing up */
    {
        if (arf_cmp_d(arb_midref(acb_realref(x)), -0.5) < 0)
        {
            reflect = 1;
            r = 0;
        }
        else if (arf_cmp_si(arb_midref(acb_realref(x)), 1) < 0)
        {
            reflect = 0;
            r = 1;
        }
        else
        {
            reflect = 0;
            r = 0;
        }

        n = 1;
    }
    else
    {
        acb_gamma_stirling_choose_param(&reflect, &r, &n, x, 1, 0, wp);
    }

    acb_init(t);
    acb_init(u);
    acb_init(v);

    if (reflect)
    {
        acb_sub_ui(t, x, 1, wp);
        acb_neg(t, t);
        acb_rising_ui_rec(u, t, r, wp);
        arb_const_pi(acb_realref(v), wp);
        acb_mul_arb(u, u, acb_realref(v), wp);
        acb_add_ui(t, t, r, wp);
        acb_gamma_stirling_eval(v, t, n, 0, wp);

        if (inverse)
        {
            /* rgamma(x) = gamma(1-x+r) sin(pi x) / ((rf(1-x, r) * pi) */
            acb_exp(v, v, wp);
            acb_sin_pi(t, x, wp);
            acb_mul(v, v, t, wp);
            acb_mul(y, u, v, wp);
            acb_div(y, v, u, prec);
        }
        else
        {
            /* gamma(x) = (rf(1-x, r) * pi) rgamma(1-x+r) csc(pi x) */
            acb_neg(v, v);
            acb_exp(v, v, wp);
            acb_sin_pi(t, x, wp);   /* todo: write a csc_pi function */
            acb_div(v, v, t, wp);
            acb_mul(y, v, u, prec);
        }
    }
    else
    {
        acb_add_ui(t, x, r, wp);
        acb_gamma_stirling_eval(u, t, n, 0, wp);

        if (inverse)
        {
            /* rgamma(x) = rf(x,r) rgamma(x+r) */
            acb_neg(u, u);
            acb_exp(u, u, prec);
            acb_rising_ui_rec(v, x, r, wp);
            acb_mul(y, v, u, prec);
        }
        else
        {
            /* gamma(x) = gamma(x+r) / rf(x,r) */
            acb_exp(u, u, prec);
            acb_rising_ui_rec(v, x, r, wp);
            acb_div(y, u, v, prec);
        }
    }

    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
}

void
acb_gamma(acb_t y, const acb_t x, slong prec)
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
acb_rgamma(acb_t y, const acb_t x, slong prec)
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
        const acb_t t_wrong, const acb_t z, ulong r, slong prec)
{
    acb_t f;
    arb_t pi, u, v;
    fmpz_t pi_mult;
    slong i, argprec;

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
acb_lgamma(acb_t y, const acb_t x, slong prec)
{
    int reflect;
    slong r, n, wp;
    acb_t t, u, v;

    if (acb_is_real(x) && arb_is_positive(acb_realref(x)))
    {
        arb_lgamma(acb_realref(y), acb_realref(x), prec);
        arb_zero(acb_imagref(y));
        return;
    }

    wp = prec + FLINT_BIT_COUNT(prec);

    acb_gamma_stirling_choose_param(&reflect, &r, &n, x, 1, 0, wp);

    acb_init(t);
    acb_init(u);
    acb_init(v);

    if (reflect)
    {
        /* log gamma(x) = log rf(1-x, r) - log gamma(1-x+r) - log sin(pi x) + log(pi) */
        acb_sub_ui(u, x, 1, wp);
        acb_neg(u, u);

        acb_rising_ui_rec(t, u, r, prec);
        acb_log(t, t, wp);
        _acb_log_rising_correct_branch(t, t, u, r, wp);

        acb_add_ui(u, u, r, wp);
        acb_gamma_stirling_eval(v, u, n, 0, wp);
        acb_sub(t, t, v, wp);

        acb_log_sin_pi(u, x, wp);
        acb_sub(t, t, u, wp);

        acb_const_pi(u, wp);
        acb_log(u, u, wp);

        acb_add(y, t, u, wp);
    }
    else
    {
        /* log gamma(x) = log gamma(x+r) - log rf(x,r) */

        acb_add_ui(t, x, r, wp);
        acb_gamma_stirling_eval(u, t, n, 0, wp);

        acb_rising_ui_rec(t, x, r, prec);
        acb_log(t, t, wp);
        _acb_log_rising_correct_branch(t, t, x, r, wp);

        acb_sub(y, u, t, prec);
    }

    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
}

