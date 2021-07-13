/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "acb.h"
#include "bernoulli.h"
#include "hypgeom.h"

/* tuning factor */
#define GAMMA_STIRLING_BETA 0.27

#define PI 3.1415926535897932385

static slong
choose_n(double log2z, double argz, int digamma, slong prec)
{
    double argf, boundn;
    slong n;

    argf = 1.0 / cos(0.5 * argz);
    argf = log(argf) * (1. / log(2));

    for (n = 1; ; n++)
    {
        if (digamma)
            boundn = bernoulli_bound_2exp_si(2*n) - (2*n)*log2z + (2*n+1)*argf;
        else
            boundn = bernoulli_bound_2exp_si(2*n) - (2*n-1)*log2z + (2*n)*argf;

        /* success */
        if (boundn <= -prec)
            return n;

        /* if the term magnitude does not decrease, r is too small */
        if (boundn > 1)
        {
            flint_printf("exception: gamma_stirling_choose_param failed to converge\n");
            flint_abort();
        }
    }
}

static void
choose_small(int * reflect, slong * r, slong * n,
    double x, double y, int use_reflect, int digamma, slong prec)
{
    double w, argz, log2z;
    slong rr;

    /* use reflection formula if very negative */
    if (x < -5.0 && use_reflect)
    {
        *reflect = 1;
        x = 1.0 - x;
    }
    else
    {
        *reflect = 0;
    }

    /* argument reduction until |z| >= w */
    w = FLINT_MAX(1.0, GAMMA_STIRLING_BETA * prec);

    rr = 0;
    while (x < 1.0 || x*x + y*y < w*w)
    {
        x++;
        rr++;
    }

    log2z = 0.5 * log(x*x + y*y) * 1.44269504088896341;
    argz = atan2(y, x);

    *r = rr;
    *n = choose_n(log2z, argz, digamma, prec);
}

static void
choose_large(int * reflect, slong * r, slong * n,
    const arf_t a, const arf_t b, int use_reflect, int digamma, slong prec)
{
    if (use_reflect && arf_sgn(a) < 0)
        *reflect = 1;
    else
        *reflect = 0;

    *r = 0;

    /* so big that we will certainly have n = 0 */
    if (arf_cmpabs_2exp_si(a, WORD_MAX / 8) >= 0 ||
        arf_cmpabs_2exp_si(b, WORD_MAX / 8) >= 0)
    {
        *n = 0;
    }
    else
    {
        slong ab, bb;
        double log2z, argz;

        ab = arf_abs_bound_lt_2exp_si(a);
        bb = arf_abs_bound_lt_2exp_si(b);

        log2z = FLINT_MAX(ab, bb);

        /* piecewise approximation of the argument */
        if (arf_is_zero(b))
        {
            if ((arf_sgn(a) < 0) && !(*reflect))
                argz = PI;
            else
                argz = 0.0;
        }
        else
        {
            if ((arf_sgn(a) < 0) && !(*reflect))
                if (arf_cmpabs(a, b) <= 0)
                    argz = PI * 0.75;
                else
                    argz = PI;
            else
                if (arf_cmpabs(a, b) <= 0)
                    argz = PI * 0.25;
                else
                    argz = PI * 0.5;
        }

        if (argz == PI)
            *n = 0;
        else
            *n = choose_n(log2z, argz, digamma, prec);
    }
}


void
acb_gamma_stirling_choose_param(int * reflect, slong * r, slong * n,
    const acb_t z, int use_reflect, int digamma, slong prec)
{
    const arf_struct * a = arb_midref(acb_realref(z));
    const arf_struct * b = arb_midref(acb_imagref(z));

    if (!arf_is_finite(a) || !arf_is_finite(b))
    {
        *reflect = *r = *n = 0;
    }
    else if (arf_cmpabs_2exp_si(a, 40) > 0 || arf_cmpabs_2exp_si(b, 40) > 0)
    {
        choose_large(reflect, r, n, a, b, use_reflect, digamma, prec);
    }
    else
    {
        choose_small(reflect, r, n,
            arf_get_d(a, ARF_RND_UP),
            arf_get_d(b, ARF_RND_UP), use_reflect, digamma, prec);
    }
}

void
arb_gamma_stirling_choose_param(int * reflect, slong * r, slong * n,
    const arb_t x, int use_reflect, int digamma, slong prec)
{
    const arf_struct * a = arb_midref(x);

    if (arf_is_inf(a) || arf_is_nan(a))
    {
        *reflect = *r = *n = 0;
    }
    else if (arf_cmpabs_2exp_si(a, 40) > 0)
    {
        arf_t b;
        arf_init(b);
        choose_large(reflect, r, n, a, b, use_reflect, digamma, prec);
        arf_clear(b);
    }
    else
    {
        choose_small(reflect, r, n,
            arf_get_d(a, ARF_RND_UP), 0.0, use_reflect, digamma, prec);
    }
}

void
acb_gamma_bound_phase(mag_t bound, const acb_t z)
{
    arf_t x, y, t, u;
    int xsign;
    slong prec;

    arf_init(x);
    arf_init(y);
    arf_init(t);
    arf_init(u);

    prec = MAG_BITS;

    /* first compute x, y such that |arg(z)| <= arg(x+yi) */

    /* argument increases with smaller real parts */
    arf_set_mag(x, arb_radref(acb_realref(z)));
    arf_sub(x, arb_midref(acb_realref(z)), x, prec, ARF_RND_FLOOR);

    xsign = arf_sgn(x);

    if (xsign >= 0)  /* argument increases away from the real axis */
        arb_get_abs_ubound_arf(y, acb_imagref(z), prec);
    else  /* argument increases closer to the real axis */
        arb_get_abs_lbound_arf(y, acb_imagref(z), prec);

    if (arf_is_zero(y))
    {
        if (xsign > 0)
            mag_one(bound);
        else
            mag_inf(bound);
    }
    else
    {
        if (xsign >= 0)
        {
            /* compute upper bound for t = y / (sqrt(x^2 + y^2) + x) */
            arf_mul(t, x, x, prec, ARF_RND_DOWN);
            arf_mul(u, y, y, prec, ARF_RND_DOWN);
            arf_add(t, t, u, prec, ARF_RND_DOWN);
            arf_sqrt(t, t, prec, ARF_RND_DOWN);
            arf_add(t, t, x, prec, ARF_RND_DOWN);
            arf_div(t, y, t, prec, ARF_RND_UP);
        }
        else
        {
            /* compute upper bound for t = (sqrt(x^2 + y^2) - x) / y */
            arf_mul(t, x, x, prec, ARF_RND_UP);
            arf_mul(u, y, y, prec, ARF_RND_UP);
            arf_add(t, t, u, prec, ARF_RND_UP);
            arf_sqrt(t, t, prec, ARF_RND_UP);
            arf_sub(t, t, x, prec, ARF_RND_UP);
            arf_div(t, t, y, prec, ARF_RND_UP);
        }

        /* compute upper bound for sqrt(1 + t^2) */
        arf_mul(t, t, t, prec, ARF_RND_UP);
        arf_add_ui(t, t, 1, prec, ARF_RND_UP);
        arf_sqrt(t, t, prec, ARF_RND_UP);

        arf_get_mag(bound, t);
    }

    arf_clear(x);
    arf_clear(y);
    arf_clear(t);
    arf_clear(u);
}

/*
  2 |B_{2n}| G(2n+k-1) / (G(k+1) G(2n+1)) |z| (T |z|^{-1})^(2n+k)
  TODO: CHECK n >= 1 ?
*/
void
acb_gamma_stirling_bound(mag_ptr err, const acb_t z, slong k0, slong knum, slong n)
{
    mag_t c, t, u, v;
    slong i, k;

    if (arb_contains_zero(acb_imagref(z)) &&
        arb_contains_nonpositive(acb_realref(z)))
    {
        for (i = 0; i < knum; i++)
            mag_inf(err + i);
        return;
    }

    mag_init(c);
    mag_init(t);
    mag_init(u);
    mag_init(v);

    /* t = lower bound for |z| */
    acb_get_mag_lower(t, z);
    /* v = upper bound for |z| */
    acb_get_mag(v, z);

    /* c = upper bound for 1/(cos(arg(z)/2) |z|) */
    acb_gamma_bound_phase(c, z);
    mag_div(c, c, t);

    /* numerator: 2 B_{2n} gamma(2n+k-1) |z| */
    mag_bernoulli_div_fac_ui(err, 2 * n);
    mag_mul_2exp_si(err, err, 1);
    mag_fac_ui(u, 2 * n + k0 - 2);
    mag_mul(err, err, u);
    mag_mul(err, err, v);

    /* denominator gamma(k+1) gamma(2n+1) */
    mag_rfac_ui(t, k0);
    mag_mul(err, err, t);

    /* multiply by c^(2n+k) */
    mag_pow_ui(t, c, 2 * n + k0);
    mag_mul(err, err, t);

    for (i = 1; i < knum; i++)
    {
        /* recurrence factor: c * (2n+k-2) / k */
        k = k0 + i;
        mag_mul(err + i, err + i - 1, c);
        mag_mul_ui(err + i, err + i, 2 * n + k - 2);
        mag_div_ui(err + i, err + i, k);
    }

    mag_clear(c);
    mag_clear(t);
    mag_clear(u);
    mag_clear(v);
}

void
arb_gamma_stirling_bound(mag_ptr err, const arb_t x, slong k0, slong knum, slong n)
{
    acb_t z;
    acb_init(z);
    acb_set_arb(z, x);
    acb_gamma_stirling_bound(err, z, k0, knum, n);
    acb_clear(z);
}

void
arb_gamma_stirling_coeff(arb_t b, ulong k, int digamma, slong prec)
{
    fmpz_t d;
    fmpz_init(d);

    BERNOULLI_ENSURE_CACHED(2 * k);

    arb_set_round_fmpz(b, fmpq_numref(bernoulli_cache + 2 * k), prec);

    if (digamma)
        fmpz_mul_ui(d, fmpq_denref(bernoulli_cache + 2 * k), 2 * k);
    else
        fmpz_mul2_uiui(d, fmpq_denref(bernoulli_cache + 2 * k), 2 * k, 2 * k - 1);

    arb_div_fmpz(b, b, d, prec);
    fmpz_clear(d);
}

void
arb_gamma_stirling_eval(arb_t s, const arb_t z, slong nterms, int digamma, slong prec)
{
    arb_t b, t, logz, zinv, zinv2;
    mag_t err;

    slong k, term_prec;
    double z_mag, term_mag;

    arb_init(b);
    arb_init(t);
    arb_init(logz);
    arb_init(zinv);
    arb_init(zinv2);

    arb_log(logz, z, prec);
    arb_inv(zinv, z, prec);

    nterms = FLINT_MAX(nterms, 1);

    arb_zero(s);

    if (nterms > 1)
    {
        arb_mul(zinv2, zinv, zinv, prec);

        z_mag = arf_get_d(arb_midref(logz), ARF_RND_UP) * 1.44269504088896;

        for (k = nterms - 1; k >= 1; k--)
        {
            term_mag = bernoulli_bound_2exp_si(2 * k);
            term_mag -= (2 * k - 1) * z_mag;
            term_prec = prec + term_mag;
            term_prec = FLINT_MIN(term_prec, prec);
            term_prec = FLINT_MAX(term_prec, 10);

            if (prec > 2000)
            {
                arb_set_round(t, zinv2, term_prec);
                arb_mul(s, s, t, term_prec);
            }
            else
                arb_mul(s, s, zinv2, term_prec);

            arb_gamma_stirling_coeff(b, k, digamma, term_prec);
            arb_add(s, s, b, term_prec);
        }

        if (digamma)
            arb_mul(s, s, zinv2, prec);
        else
            arb_mul(s, s, zinv, prec);
    }

    /* remainder bound */
    mag_init(err);
    arb_gamma_stirling_bound(err, z, digamma ? 1 : 0, 1, nterms);
    mag_add(arb_radref(s), arb_radref(s), err);
    mag_clear(err);

    if (digamma)
    {
        arb_neg(s, s);
        arb_mul_2exp_si(zinv, zinv, -1);
        arb_sub(s, s, zinv, prec);
        arb_add(s, s, logz, prec);
    }
    else
    {
        /* (z-0.5)*log(z) - z + log(2*pi)/2 */
        arb_one(t);
        arb_mul_2exp_si(t, t, -1);
        arb_sub(t, z, t, prec);
        arb_mul(t, logz, t, prec);
        arb_add(s, s, t, prec);
        arb_sub(s, s, z, prec);
        arb_const_log_sqrt2pi(t, prec);
        arb_add(s, s, t, prec);
    }

    arb_clear(t);
    arb_clear(b);
    arb_clear(zinv);
    arb_clear(zinv2);
    arb_clear(logz);
}

void
arb_gamma_fmpq_stirling(arb_t y, const fmpq_t a, slong prec)
{
    int reflect;
    slong r, n, wp;
    arb_t x, t, u, v;

    wp = prec + FLINT_BIT_COUNT(prec);

    arb_init(x);
    arb_init(t);
    arb_init(u);
    arb_init(v);

    arb_set_fmpq(x, a, wp);
    arb_gamma_stirling_choose_param(&reflect, &r, &n, x, 1, 0, wp);

    if (reflect)
    {
        /* gamma(x) = (rf(1-x, r) * pi) / (gamma(1-x+r) sin(pi x)) */
        fmpq_t b;
        fmpq_init(b);
        fmpz_sub(fmpq_numref(b), fmpq_denref(a), fmpq_numref(a));
        fmpz_set(fmpq_denref(b), fmpq_denref(a));
        arb_rising_fmpq_ui(u, b, r, wp);
        fmpq_clear(b);
        arb_const_pi(v, wp);
        arb_mul(u, u, v, wp);
        arb_sub_ui(t, x, 1, wp);
        arb_neg(t, t);
        arb_add_ui(t, t, r, wp);
        arb_gamma_stirling_eval(v, t, n, 0, wp);
        arb_exp(v, v, wp);
        arb_sin_pi_fmpq(t, a, wp);
        arb_mul(v, v, t, wp);
    }
    else
    {
        /* gamma(x) = gamma(x+r) / rf(x,r) */
        arb_add_ui(t, x, r, wp);
        arb_gamma_stirling_eval(u, t, n, 0, wp);
        arb_exp(u, u, prec);
        arb_rising_fmpq_ui(v, a, r, wp);
    }

    arb_div(y, u, v, prec);

    arb_clear(t);
    arb_clear(u);
    arb_clear(v);
    arb_clear(x);
}

void
arb_gamma_const_1_3_eval(arb_t s, slong prec)
{
    hypgeom_t series;
    arb_t t, u;
    slong wp = prec + 4 + 2 * FLINT_BIT_COUNT(prec);

    arb_init(t);
    arb_init(u);

    hypgeom_init(series);

    fmpz_poly_set_str(series->A, "1  1");
    fmpz_poly_set_str(series->B, "1  1");
    fmpz_poly_set_str(series->P, "4  5 -46 108 -72");
    fmpz_poly_set_str(series->Q, "4  0 0 0 512000");

    prec += FLINT_CLOG2(prec);

    arb_hypgeom_infsum(s, t, series, wp, wp);

    arb_sqrt_ui(u, 10, wp);
    arb_mul(t, t, u, wp);

    arb_const_pi(u, wp);
    arb_pow_ui(u, u, 4, wp);
    arb_mul_ui(u, u, 12, wp);
    arb_mul(s, s, u, wp);

    arb_div(s, s, t, wp);
    arb_root_ui(s, s, 2, wp);
    arb_root_ui(s, s, 3, prec);

    hypgeom_clear(series);
    arb_clear(t);
    arb_clear(u);
}

ARB_DEF_CACHED_CONSTANT(arb_gamma_const_1_3, arb_gamma_const_1_3_eval)

void
arb_gamma_const_1_4_eval(arb_t x, slong prec)
{
    arb_t t, u;
    slong wp = prec + 4 + 2 * FLINT_BIT_COUNT(prec);

    arb_init(t);
    arb_init(u);

    arb_one(t);
    arb_sqrt_ui(u, 2, wp);
    arb_agm(x, t, u, wp);

    arb_const_pi(t, wp);
    arb_mul_2exp_si(t, t, 1);
    arb_sqrt(u, t, wp);
    arb_mul(u, u, t, wp);

    arb_div(x, u, x, wp);
    arb_sqrt(x, x, wp);

    arb_clear(t);
    arb_clear(u);
}

ARB_DEF_CACHED_CONSTANT(arb_gamma_const_1_4, arb_gamma_const_1_4_eval)

void
arb_gamma_small_frac(arb_t y, unsigned int p, unsigned int q, slong prec)
{
    slong wp = prec + 4;

    if (q == 1)
    {
        arb_one(y);
    }
    else if (q == 2)  /* p = 1 */
    {
        arb_const_sqrt_pi(y, prec);
    }
    else if (q == 3)
    {
        if (p == 1)
        {
            arb_gamma_const_1_3(y, prec);
        }
        else  /* p = 2 */
        {
            arb_t t;
            arb_init(t);
            arb_gamma_const_1_3(y, wp);
            arb_sqrt_ui(t, 3, wp);
            arb_mul(y, y, t, wp);
            arb_const_pi(t, wp);
            arb_div(y, t, y, prec);
            arb_mul_2exp_si(y, y, 1);
            arb_clear(t);
        }
    }
    else if (q == 4)
    {
        if (p == 1)
        {
            arb_gamma_const_1_4(y, prec);
        }
        else  /* p = 3 */
        {
            arb_t t;
            arb_init(t);
            arb_sqrt_ui(y, 2, wp);
            arb_const_pi(t, wp);
            arb_mul(y, y, t, wp);
            arb_gamma_const_1_4(t, wp);
            arb_div(y, y, t, prec);
            arb_clear(t);
        }
    }
    else if (q == 6)
    {
        arb_t t;
        arb_init(t);
        arb_const_pi(t, wp);
        arb_div_ui(t, t, 3, wp);
        arb_sqrt(t, t, wp);
        arb_set_ui(y, 2);
        arb_root_ui(y, y, 3, wp);
        arb_mul(t, t, y, wp);
        arb_gamma_const_1_3(y, wp);
        arb_mul(y, y, y, prec);

        if (p == 1)
        {
            arb_div(y, y, t, prec);
        }
        else  /* p = 5 */
        {
            arb_div(y, t, y, wp);
            arb_const_pi(t, wp);
            arb_mul(y, y, t, prec);
            arb_mul_2exp_si(y, y, 1);
        }

        arb_clear(t);
    }
    else
    {
        flint_printf("small fraction not implemented!\n");
        flint_abort();
    }
}

slong _arb_compute_bs_exponents(slong * tab, slong n);
slong _arb_get_exp_pos(const slong * tab, slong step);

static void
bsplit2(arb_t P, arb_t Q, const fmpz_t zp, const fmpz_t zq,
    const slong * xexp, arb_srcptr xpow,
    ulong N, slong a, slong b, int cont, slong prec)
{
    if (b - a == 1)
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_set(t, zp);
        fmpz_addmul_ui(t, zq, a + 1);
        arb_set_fmpz(P, t);
        arb_set(Q, P);
        fmpz_clear(t);
    }
    else
    {
        arb_t Pb, Qb;
        slong step, i, m;

        arb_init(Pb);
        arb_init(Qb);

        step = (b - a) / 2;
        m = a + step;

        bsplit2(P, Q, zp, zq, xexp, xpow, N, a, m, 1, prec);
        bsplit2(Pb, Qb, zp, zq, xexp, xpow, N, m, b, 1, prec);

        arb_mul(P, P, Pb, prec);
        arb_mul(Q, Q, Pb, prec);

        i = (step == 1) ? 0 : _arb_get_exp_pos(xexp, step);
        arb_addmul(Q, Qb, xpow + i, prec);

        arb_clear(Pb);
        arb_clear(Qb);
    }
}

static void
bsplit3(arb_t P, arb_t Q, const fmpz_t zp, const fmpz_t zq,
    const slong * xexp, arb_srcptr xpow,
    ulong N, slong a, slong b, int cont, slong prec)
{
    if (b - a == 1)
    {
        fmpz_t t;
        fmpz_init(t);
        arb_set(P, xpow + 0);  /* N zq */
        fmpz_set(t, zp);
        fmpz_submul_ui(t, zq, a + 1);  /* zp - (a + 1) zq */
        arb_set_fmpz(Q, t);
        fmpz_clear(t);
    }
    else
    {
        arb_t Pb, Qb;
        slong step, i, m;

        arb_init(Pb);
        arb_init(Qb);

        step = (b - a) / 2;
        m = a + step;

        bsplit3(P, Q, zp, zq, xexp, xpow, N, a, m, 1, prec);
        bsplit3(Pb, Qb, zp, zq, xexp, xpow, N, m, b, 1, prec);

        i = _arb_get_exp_pos(xexp, m - a);

        arb_mul(P, P, xpow + i, prec);
        if (b - m != m - a)
            arb_mul(P, P, xpow + 0, prec);

        arb_addmul(P, Q, Pb, prec);

        if (cont)
        {
            arb_mul(Q, Q, Qb, prec);
        }
        else
        {
            i = _arb_get_exp_pos(xexp, m - a);

            arb_mul(Q, xpow + i, xpow + i, prec);
            if (b - m != m - a)
                arb_mul(Q, Q, xpow + 0, prec);
        }

        arb_clear(Pb);
        arb_clear(Qb);
    }
}

double d_lambertw_branch1(double x);

static ulong
more_trailing_zeros(ulong N)
{
    ulong bc, N2;

    bc = FLINT_BIT_COUNT(N);

    if (bc >= 8)
    {
        N2 = (N >> (bc - 5)) << (bc - 5);
        N2 += ((N2 != N) << (bc - 5));
        return N2;
    }
    else
    {
        return N;
    }
}

#define C_LOG2 0.69314718055994530942
#define C_EXP1 2.7182818284590452354

static void
build_bsplit_power_table(arb_ptr xpow, const slong * xexp, slong len, slong prec)
{
    slong i;

    for (i = 1; i < len; i++)
    {
        if (xexp[i] == 2 * xexp[i-1])
        {
            arb_mul(xpow + i, xpow + i - 1, xpow + i - 1, prec);
        }
        else if (xexp[i] == 2 * xexp[i-2]) /* prefer squaring if possible */
        {
            arb_mul(xpow + i, xpow + i - 2, xpow + i - 2, prec);
        }
        else if (xexp[i] == 2 * xexp[i-1] + 1)
        {
            arb_mul(xpow + i, xpow + i - 1, xpow + i - 1, prec);
            arb_mul(xpow + i, xpow + i, xpow, prec);
        }
        else if (xexp[i] == 2 * xexp[i-2] + 1)
        {
            arb_mul(xpow + i, xpow + i - 2, xpow + i - 2, prec);
            arb_mul(xpow + i, xpow + i, xpow, prec);
        }
        else
        {
            flint_printf("power table has the wrong structure!\n");
            flint_abort();
        }
    }
}

/* assumes z in [1, 2] */
static void
arb_gamma_fmpq_general_off1(arb_t res, const fmpq_t z, slong prec)
{
    slong wp, N, n, n2, length, length2, wp2;
    double alpha;
    arb_t P, Q;
    slong *xexp, *xexp2;
    arb_ptr xpow;
    mag_t err, err2;

    wp = prec + 30;

    alpha = 0.52;  /* tuning parameter between 0.5 and 1 */

    N = alpha * C_LOG2 * wp;
    N = more_trailing_zeros(N);
    alpha = N / (C_LOG2 * wp);

    /* Terms in convergent series */
    n = (1 - alpha) / d_lambertw((1 - alpha) / (alpha * C_EXP1)) * C_LOG2 * wp;

    /* Precision and terms in asymptotic series */
    wp2 = wp * (1 - alpha);
    wp2 = FLINT_MAX(wp2, 30);
    n2 = (alpha - 1) / d_lambertw_branch1((alpha - 1) / (alpha * C_EXP1)) * C_LOG2 * wp;
    n2 = FLINT_MAX(n2, 2);  /* binary splitting correctness */

    mag_init(err);
    mag_init(err2);
    arb_init(P);
    arb_init(Q);

    /* compute the powers of x = N*zq that will appear (at least x^1) */
    xexp = flint_calloc(2 * FLINT_BITS, sizeof(slong));
    xexp2 = flint_calloc(2 * FLINT_BITS, sizeof(slong));

    length = _arb_compute_bs_exponents(xexp, n);
    length2 = _arb_compute_bs_exponents(xexp2, n2);

    xpow = _arb_vec_init(FLINT_MAX(length, length2));

    arb_set_fmpz(xpow + 0, fmpq_denref(z));
    arb_mul_ui(xpow + 0, xpow + 0, N, wp);

    build_bsplit_power_table(xpow, xexp, length, wp);

    /* 1F1(1, 1+z, N) */
    bsplit2(P, Q, fmpq_numref(z), fmpq_denref(z), xexp, xpow, N, 0, n, 0, wp);
    arb_div(P, Q, P, wp);

    /* Convergent series error bound: N^n / n! (1 + (N/n) + ...) */
    mag_set_ui(err, N);
    mag_pow_ui(err, err, n);
    mag_rfac_ui(err2, n);
    mag_mul(err, err, err2);
    mag_set_ui(err2, N);
    mag_div_ui(err2, err2, n);
    mag_geom_series(err2, err2, 0);
    mag_mul(err, err, err2);
    arb_add_error_mag(P, err);

    /* divide 1F1 by z */
    arb_mul_fmpz(P, P, fmpq_denref(z), wp);
    arb_div_fmpz(P, P, fmpq_numref(z), wp);
    arb_swap(res, P);

    build_bsplit_power_table(xpow, xexp2, length2, wp2);

    bsplit3(P, Q, fmpq_numref(z), fmpq_denref(z), xexp2, xpow, N, 0, n2, 0, wp2);
    arb_div(P, P, Q, wp2);

    /* 2F0 error bound (bounded by first omitted term) */
    mag_fac_ui(err, n2);
    mag_set_ui_lower(err2, N);
    mag_pow_ui_lower(err2, err2, n2);
    mag_div(err, err, err2);
    arb_add_error_mag(P, err);

    /* N^z * exp(-N) * (1F1/z + 2F0/N)  */
    arb_div_ui(P, P, N, wp2);

    arb_add(res, res, P, wp);
    arb_set_ui(Q, N);
    arb_pow_fmpq(Q, Q, z, wp);
    arb_mul(res, res, Q, wp);

    arb_set_si(Q, -N);
    arb_exp(Q, Q, wp);
    arb_mul(res, res, Q, wp);

    _arb_vec_clear(xpow, FLINT_MAX(length, length2));
    flint_free(xexp);
    flint_free(xexp2);

    arb_clear(P);
    arb_clear(Q);
    mag_clear(err);
    mag_clear(err2);
}

/* assumes z in (0, 1] */
void
arb_gamma_fmpq_hyp(arb_t res, const fmpq_t z, slong prec)
{
    fmpq_t t;
    fmpq_init(t);
    fmpq_add_ui(t, z, 1);
    arb_gamma_fmpq_general_off1(res, t, prec);
    arb_mul_fmpz(res, res, fmpq_denref(z), prec + 30);
    arb_div_fmpz(res, res, fmpq_numref(z), prec);
    fmpq_clear(t);
}

void
arb_gamma_fmpq_outward(arb_t y, const fmpq_t x, slong prec)
{
    fmpq_t a;
    fmpz_t n;
    fmpz p, q;
    slong m;
    arb_t t, u;

    fmpq_init(a);
    fmpz_init(n);
    arb_init(t);
    arb_init(u);

    /* write x = a + n with 0 < a <= 1 */
    if (fmpz_is_one(fmpq_denref(x)))
    {
        fmpq_one(a);
        fmpz_sub_ui(n, fmpq_numref(x), 1);
    }
    else
    {
        fmpz_fdiv_qr(n, fmpq_numref(a), fmpq_numref(x), fmpq_denref(x));
        fmpz_set(fmpq_denref(a), fmpq_denref(x));
    }

    if (!fmpz_fits_si(n))
    {
        flint_printf("gamma: too large fmpq to reduce to 0!\n");
        flint_abort();
    }

    m = fmpz_get_si(n);

    /* evaluate gamma(a) */
    p = *fmpq_numref(a);
    q = *fmpq_denref(a);

    if (q == 1 || q == 2 || q == 3 || q == 4 || q == 6)
    {
        arb_gamma_small_frac(t, p, q, prec);
    }
    else
    {
        arb_gamma_fmpq_hyp(t, a, prec);
    }

    /* argument reduction */
    if (m >= 0)
    {
        arb_rising_fmpq_ui(u, a, m, prec);
        arb_mul(y, t, u, prec);
    }
    else
    {
        arb_rising_fmpq_ui(u, x, -m, prec);
        arb_div(y, t, u, prec);
    }

    fmpq_clear(a);
    fmpz_clear(n);
    arb_clear(t);
    arb_clear(u);
}

void
arb_gamma_fmpq(arb_t y, const fmpq_t x, slong prec)
{
    fmpz p, q;

    p = *fmpq_numref(x);
    q = *fmpq_denref(x);

    if ((q == 1 || q == 2 || q == 3 || q == 4 || q == 6) && !COEFF_IS_MPZ(p))
    {
        if (q == 1)
        {
            if (p <= 0)
            {
                arb_indeterminate(y);
                return;
            }

            if (p < 1200 || 1.44265 * (p*log(p) - p) < 15.0 * prec)
            {
                fmpz_t t;
                fmpz_init(t);
                fmpz_fac_ui(t, p - 1);
                arb_set_round_fmpz(y, t, prec);
                fmpz_clear(t);
                return;
            }
        }

        p = FLINT_ABS(p);

        if (p < q * 500.0 || p < q * (500.0 + 0.1 * prec * sqrt(prec)))
        {
            arb_gamma_fmpq_outward(y, x, prec);
            return;
        }
    }

    if (q != 1 && prec > 7000 + 300 * fmpz_bits(fmpq_denref(x)) && 
        (slong) fmpz_bits(fmpq_numref(x)) - (slong) fmpz_bits(fmpq_denref(x)) < FLINT_BITS &&
        fabs(fmpq_get_d(x)) < 0.03 * prec * sqrt(prec))
    {
        arb_gamma_fmpq_outward(y, x, prec);
        return;
    }

    arb_gamma_fmpq_stirling(y, x, prec);
}

void
arb_gamma_fmpz(arb_t y, const fmpz_t x, slong prec)
{
    fmpq_t t;
    *fmpq_numref(t) = *x;
    *fmpq_denref(t) = WORD(1);
    arb_gamma_fmpq(y, t, prec);
}

static void
_arb_gamma(arb_t y, const arb_t x, slong prec, int inverse)
{
    int reflect;
    slong r, n, wp;
    arb_t t, u, v;
    double acc;

    if (arb_is_exact(x))
    {
        const arf_struct * mid = arb_midref(x);

        if (arf_is_special(mid))
        {
            if (!inverse && arf_is_pos_inf(mid))
                arb_set(y, x);
            else if (arf_is_nan(mid) || arf_is_neg_inf(mid) || !inverse)
                arb_indeterminate(y);
            else
                arb_zero(y);
            return;
        }
        else if (inverse && arf_is_int(mid) && arf_sgn(mid) < 0)
        {
            arb_zero(y);
        }
        else
        {
            /* fast gamma(n), gamma(n/2) or gamma(n/4), ... */
            if (arf_cmpabs_2exp_si(mid, prec) < 0 &&
                (arf_is_int_2exp_si(mid, -2) || (prec > 1000 && arf_is_int_2exp_si(mid, -1000))))
            {
                fmpq_t a;
                fmpq_init(a);
                arf_get_fmpq(a, mid);
                arb_gamma_fmpq(y, a, prec + 2 * inverse);
                if (inverse)
                    arb_inv(y, y, prec);
                fmpq_clear(a);
                return;
            }
        }
    }

    /* todo: for large x (if exact or accurate enough), increase precision */
    acc = arb_rel_accuracy_bits(x);
    acc = FLINT_MAX(acc, 0);
    wp = FLINT_MIN(prec, acc + 20);
    wp = FLINT_MAX(wp, 2);
    wp = wp + FLINT_BIT_COUNT(wp);

    if (acc < 3)  /* try to avoid divisions blowing up */
    {
        if (arf_cmp_d(arb_midref(x), -0.5) < 0)
        {
            reflect = 1;
            r = 0;
        }
        else if (arf_cmp_si(arb_midref(x), 1) < 0)
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
        arb_gamma_stirling_choose_param(&reflect, &r, &n, x, 1, 0, wp);
    }

    arb_init(t);
    arb_init(u);
    arb_init(v);

    if (reflect)
    {
        arb_sub_ui(t, x, 1, wp);
        arb_neg(t, t);
        arb_rising_ui_rec(u, t, r, wp);
        arb_const_pi(v, wp);
        arb_mul(u, u, v, wp);
        arb_add_ui(t, t, r, wp);
        arb_gamma_stirling_eval(v, t, n, 0, wp);

        if (inverse)
        {
            /* rgamma(x) = gamma(1-x+r) sin(pi x) / ((rf(1-x, r) * pi) */
            arb_exp(v, v, wp);
            arb_sin_pi(t, x, wp);
            arb_mul(v, v, t, wp);
            arb_mul(y, u, v, wp);
            arb_div(y, v, u, prec);
        }
        else
        {
            /* gamma(x) = (rf(1-x, r) * pi) rgamma(1-x+r) csc(pi x) */
            arb_neg(v, v);
            arb_exp(v, v, wp);
            arb_csc_pi(t, x, wp);
            arb_mul(v, v, t, wp);
            arb_mul(y, v, u, prec);
        }
    }
    else
    {
        arb_add_ui(t, x, r, wp);
        arb_gamma_stirling_eval(u, t, n, 0, wp);

        if (inverse)
        {
            /* rgamma(x) = rf(x,r) rgamma(x+r) */
            arb_neg(u, u);
            arb_exp(u, u, prec);
            arb_rising_ui_rec(v, x, r, wp);
            arb_mul(y, v, u, prec);
        }
        else
        {
            /* gamma(x) = gamma(x+r) / rf(x,r) */
            arb_exp(u, u, prec);
            arb_rising_ui_rec(v, x, r, wp);
            arb_div(y, u, v, prec);
        }
    }

    arb_clear(t);
    arb_clear(u);
    arb_clear(v);
}

void
arb_gamma(arb_t y, const arb_t x, slong prec)
{
    _arb_gamma(y, x, prec, 0);
}

void
arb_rgamma(arb_t y, const arb_t x, slong prec)
{
    _arb_gamma(y, x, prec, 1);
}

void
arb_lgamma(arb_t y, const arb_t x, slong prec)
{
    int reflect;
    slong r, n, wp;
    arb_t t, u;

    if (!arb_is_positive(x))
    {
        arb_indeterminate(y);
        return;
    }

    /* fast gamma(n), gamma(n/2) or gamma(n/4), ... */
    if (arb_is_exact(x) && arf_cmpabs_2exp_si(arb_midref(x), prec) < 0 &&
        (arf_is_int_2exp_si(arb_midref(x), -2) || (prec > 10000 && arf_is_int_2exp_si(arb_midref(x), -1000))))
    {
        fmpq_t a;
        fmpq_init(a);
        arf_get_fmpq(a, arb_midref(x));
        arb_gamma_fmpq(y, a, prec);
        arb_log(y, y, prec);
        fmpq_clear(a);
        return;
    }

    wp = prec + FLINT_BIT_COUNT(prec);

    arb_gamma_stirling_choose_param(&reflect, &r, &n, x, 0, 0, wp);

    /* log(gamma(x)) = log(gamma(x+r)) - log(rf(x,r)) */
    arb_init(t);
    arb_init(u);

    arb_add_ui(t, x, r, wp);
    arb_gamma_stirling_eval(u, t, n, 0, wp);
    arb_rising_ui_rec(t, x, r, wp);
    arb_log(t, t, wp);
    arb_sub(y, u, t, prec);

    arb_clear(t);
    arb_clear(u);
}

