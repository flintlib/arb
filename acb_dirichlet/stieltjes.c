/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "flint/double_extras.h"
#include "acb_dirichlet.h"
#include "acb_calc.h"

/* The integrand. */
int _f_stieltjes(acb_ptr res, const acb_t x, void * param, slong order, slong prec);

/* Compute the approximate saddle point. */
static void
stieltjes_omega(acb_t omega, const fmpz_t n1, slong prec)
{
    acb_t t, u;
    fmpz_t k;

    acb_init(t);
    acb_init(u);
    fmpz_init(k);

    arb_set_fmpz(acb_imagref(t), n1);
    acb_const_pi(u, prec);
    acb_mul_2exp_si(u, u, 1);
    acb_div(u, t, u, prec);
    acb_lambertw(t, u, k, 0, prec);
    acb_div(t, u, t, prec);
    acb_mul_2exp_si(t, t, 1);
    acb_sub_ui(t, t, 1, prec);
    acb_neg(t, t);
    acb_mul_onei(t, t);
    acb_mul_2exp_si(t, t, -1);

    acb_set(omega, t);

    acb_clear(t);
    acb_clear(u);
    fmpz_clear(k);
}

/* Compute an approximation of the magnitude of gamma_n (given n1 = n+1). */
static void
stieltjes_mag_approx(arb_t C, mag_t tol, const fmpz_t n1)
{
    slong prec;
    acb_t w, v, q;

    prec = 32 + 2 * fmpz_bits(n1);

    acb_init(w);
    acb_init(v);
    acb_init(q);

    stieltjes_omega(w, n1, prec);
    _f_stieltjes(v, w, (void *) n1, 0, prec);

    acb_set_fmpz(q, n1);
    acb_sqrt(q, q, prec);
    acb_mul(v, v, q, prec);

    acb_get_mag(tol, v);

    arb_set(C, acb_imagref(w));
    mag_zero(arb_radref(C));

    acb_clear(w);
    acb_clear(v);
    acb_clear(q);
}

static void
stieltjes_bound_quadratic_term(arb_t B, const acb_t z, const fmpz_t n1, slong prec)
{
    acb_t t, log_t;
    mag_t b, c, d;

    acb_init(t);
    acb_init(log_t);
    mag_init(b);
    mag_init(c);
    mag_init(d);

    /* g''(z) = (n+1)(1+1/log(t)) / (t^2 log(t)) */
    /* t = 1/2 + iz, log_t = log(t) */
    acb_mul_onei(t, z);
    acb_set_d(log_t, 0.5);
    acb_add(t, t, log_t, prec);
    acb_log(log_t, t, prec);

    acb_get_mag_lower(b, t);
    acb_get_mag_lower(c, log_t);

    /* d = 1+1/log(t) */
    mag_inv(d, c);
    mag_add_ui(d, d, 1);

    /* divide by t^2 log(t) */
    mag_div(d, d, b);
    mag_div(d, d, b);
    mag_div(d, d, c);

    mag_mul_fmpz(d, d, n1);

    /* For Taylor remainder: 1/2 rad(z)^2 */
    mag_hypot(b, arb_radref(acb_realref(z)), arb_radref(acb_imagref(z)));
    mag_mul(b, b, b);
    mag_mul_2exp_si(b, b, -1);

    mag_mul(d, d, b);

    arf_set_mag(arb_midref(B), d);
    mag_zero(arb_radref(B));

    acb_clear(t);
    acb_clear(log_t);
    mag_clear(b);
    mag_clear(c);
    mag_clear(d);
}

static void
stieltjes_bound_large3(arb_t B, const acb_t x, const fmpz_t n1, slong prec)
{
    acb_t m, t, log_t, u, v;
    acb_t g0, g1, g2;
    mag_t C, D;

    acb_init(m);
    acb_init(t);
    acb_init(log_t);
    acb_init(u);
    acb_init(v);

    acb_init(g0);
    acb_init(g1);
    acb_init(g2);

    mag_init(C);
    mag_init(D);

    /* m = mid(x) */
    acb_set(m, x);
    mag_zero(arb_radref(acb_realref(m)));
    mag_zero(arb_radref(acb_imagref(m)));

    /* t = 1/2 + im */
    acb_mul_onei(t, m);
    acb_set_d(u, 0.5);
    acb_add(t, t, u, prec);

    /* log_t = log(1/2+im) */
    acb_log(log_t, t, prec);

    /* u = x - m */
    acb_sub(u, x, m, prec);
    acb_get_mag(D, u);

    /* g0 = g(m) = (n+1) log(log(1/2+im)) - 2 pi m */
    acb_log(g0, log_t, prec);
    acb_mul_fmpz(g0, g0, n1, prec);
    acb_const_pi(v, prec);
    acb_mul(v, v, m, prec);
    acb_mul_2exp_si(v, v, 1);
    acb_sub(g0, g0, v, prec);

    /* g1 = g'(m) (x-m); g'(m) = i(n+1)/((1/2+im) log(1/2+im)) - 2 pi */
    acb_mul(g1, t, log_t, prec);
    acb_inv(g1, g1, prec);
    acb_mul_fmpz(g1, g1, n1, prec);
    acb_mul_onei(g1, g1);
    acb_const_pi(t, prec);  /* recycle t */
    acb_mul_2exp_si(t, t, 1);
    acb_sub(g1, g1, t, prec);
    acb_mul(g1, g1, u, prec);
    /* absolute value -- optional (does not seem to affect speed much) */
    acb_abs(acb_realref(g1), g1, prec);
    arb_zero(acb_imagref(g1));

    stieltjes_bound_quadratic_term(acb_realref(g2), x, n1, prec);

    acb_exp(g0, g0, prec);
    acb_exp(g1, g1, prec);
    acb_exp(g2, g2, prec);

    acb_get_mag(C, g0);
    acb_get_mag(D, g1);
    mag_mul(C, C, D);
    acb_get_mag(D, g2);
    mag_mul(C, C, D);
    mag_mul_ui(C, C, 5);

    arb_zero(B);
    arf_set_mag(arb_midref(B), C);

    acb_clear(m);
    acb_clear(t);
    acb_clear(log_t);
    acb_clear(u);
    acb_clear(v);
    mag_clear(C);
    mag_clear(D);

    acb_clear(g0);
    acb_clear(g1);
    acb_clear(g2);
}

static void
stieltjes_bound_large(acb_t res, const acb_t x, const fmpz_t n1, slong prec)
{
    arb_t B;
    mag_t t;

    arb_init(B);
    mag_init(t);

    prec = FLINT_MIN(prec, 30 + fmpz_bits(n1));

    stieltjes_bound_large3(B, x, n1, prec);
    arb_get_mag(t, B);
    acb_zero(res);
    arb_add_error_mag(acb_realref(res), t);
    arb_add_error_mag(acb_imagref(res), t);

    arb_clear(B);
    mag_clear(t);
}

int
_f_stieltjes(acb_ptr res, const acb_t x, void * param, slong order, slong prec)
{
    acb_t t, u;
    const fmpz * n1;

    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    n1 = (const fmpz *)(param);

    acb_init(t);
    acb_init(u);

    acb_mul_onei(t, x);
    acb_set_d(u, 0.5);
    acb_add(t, t, u, prec);
    acb_set_ui(u, 5);

    /* check branch cut of logarithm */
    if (arb_contains_zero(acb_imagref(t)) && !arb_is_positive(acb_realref(t)))
    {
        acb_indeterminate(res);
    }
    else if ((order == 1 || acb_rel_accuracy_bits(x) < prec - 10)
        && arb_gt(acb_realref(t), acb_realref(u))
        && arb_is_positive(acb_imagref(t)))
    {
        /* note: requires re(t) > 1, im(t) > 0 */
        stieltjes_bound_large(res, x, n1, prec);
    }
    else
    {
        acb_const_pi(u, prec);
        acb_mul(u, u, x, prec);
        acb_sech(u, u, prec);

        if (acb_is_finite(u))
        {
            acb_mul(u, u, u, prec);
            acb_log(t, t, prec);
            acb_pow_fmpz(t, t, n1, prec);
            acb_mul(res, t, u, prec);
        }
        else
        {
            acb_indeterminate(res);
        }
    }

    acb_clear(t);
    acb_clear(u);

    return 0;
}

/* Estimate log2(|gamma_n|), for fast cancellation estimate when n is small. */
static double
stieltjes_mag(double n)
{
    double u, v, A, B;
    double va, vb, t;
    double pi = 3.141592653589793;
    int i;

    va = 1e-6;
    vb = 0.5 * pi - 1e-6;

    for (i = 0; i < 53; i++)  /* bisection */
    {
        v = (va + vb) * 0.5;
        t = 2 * pi * exp(v * tan(v)) - n * cos(v) / v;
        if (t < 0.0)
            va = v;
        else
            vb = v;
    }

    v = va;
    u = v * tan(v);
    A = 0.5 * log(u * u + v * v) - u / (u * u + v * v);
    B = 2 * sqrt(2 * pi) * sqrt(u * u + v * v) * pow((u + 1) * (u + 1) + v * v, -0.25);
    t = (log(B) + n * A - 0.5 * log(n)) / log(2);
    return t;
}

static double _hypot(double x, double y) { return sqrt(x * x + y * y); }

/* log2 magnitude of integrand at x+yi */
static double
integrand_mag(double n, double x, double y)
{
    double t, u;
    t = log(_hypot(0.5 - y, x));
    u = atan2(x, 0.5 - y);
    t = log(_hypot(t,u)) * (n+1) - 2.0 * 3.1415926535897932 * x;
    return t * 1.44269504088896341;
}

static double
find_x_maximizing_mag(double n, double y)
{
    double xa, xb, xma, xmb, ma, mb;
    int i;

    xa = 1.0;
    xb = n;

    for (i = 0; i < 80; i++)
    {
        xma = xa + (xb - xa) / 3;
        xmb = xb - (xb - xa) / 3;

        ma = integrand_mag(n, xma, y);
        mb = integrand_mag(n, xmb, y);

        if (ma < mb)
            xa = xma;
        else
            xb = xmb;
    }

    return xa;
}

static void
stieltjes_choose_N(arb_t N, const fmpz_t n1, slong prec)
{
    if (fmpz_bits(n1) < 30)
    {
        double nn, NN;
        nn = fmpz_get_d(n1) - 1.0;
        NN = FLINT_MAX(nn, 4);
        while (integrand_mag(nn, NN, 0.0) > -prec - 20)
            NN *= 2.0;
        arb_set_d(N, NN);
    }
    else
    {
        arb_set_fmpz(N, n1);
    }
}

static void
stieltjes_tail_bound(mag_t bound, const arb_t N, const fmpz_t n1)
{
    slong prec;
    arb_t x, y;
    arb_init(x);
    arb_init(y);
    prec = MAG_BITS + fmpz_bits(n1);
    arb_set(x, N);
    arb_mul_2exp_si(x, x, 1);
    arb_const_pi(y, prec);
    arb_mul(y, y, x, prec);
    arb_neg(y, y);
    arb_exp(y, y, prec);
    arb_log(x, x, prec);
    arb_pow_fmpz(x, x, n1, prec);
    arb_mul(x, x, y, prec);
    arb_get_mag(bound, x);
    arb_clear(x);
    arb_clear(y);
}

void
acb_dirichlet_stieltjes_integral(arb_t res, const fmpz_t n1, slong prec)
{
    double gamma_mag, max_mag, cancellation, xa;
    acb_t a, b, v, w;
    slong wp;
    mag_t tol, bound;
    acb_calc_integrate_opt_t opt;
    /* integration points */
    arb_t M, N, C;

    arb_init(M);
    arb_init(N);
    arb_init(C);

    acb_init(a);
    acb_init(b);
    acb_init(v);
    acb_init(w);
    mag_init(tol);
    mag_init(bound);

    arb_set_ui(M, 10);

    stieltjes_choose_N(N, n1, prec);
    stieltjes_tail_bound(bound, N, n1);

    if (fmpz_cmp_ui(n1, 5000) < 0)
    {
        double nn = fmpz_get_ui(n1) - 1;
        gamma_mag = stieltjes_mag(nn);
        xa = find_x_maximizing_mag(nn, 0.0);
        max_mag = integrand_mag(nn, xa, 0.0);
        cancellation = FLINT_MAX(0, max_mag - gamma_mag);

        if (cancellation < 10 + 0.1 * prec)
        {
            arb_zero(C);
            mag_one(tol);
            mag_mul_2exp_si(tol, tol, gamma_mag);
        }
        else
        {
            stieltjes_mag_approx(C, tol, n1);
            cancellation = 0;
        }
    }
    else
    {
        stieltjes_mag_approx(C, tol, n1);
        cancellation = 0;
    }

    mag_mul_2exp_si(tol, tol, -prec - 5);

    /* todo: 1 * fmpz_bits(n1) should be enough, but acb powering is
       inaccurate */
    wp = prec + 2 * fmpz_bits(n1) + cancellation + 10;

    acb_calc_integrate_opt_init(opt);

    /* opt->verbose = 1; */
    opt->deg_limit = 100 + 1.2 * prec;   /* Small speedup. */

    if (arb_is_zero(C))
    {
        acb_zero(a);
        acb_set_arb(b, N);
        acb_calc_integrate(w, _f_stieltjes, (void *) n1, a, b, wp, tol, opt, wp);
        acb_add(v, v, w, wp);
    }
    else
    {
        acb_zero(a);                 /* a = 0 */
        acb_set_arb(b, M);           /* b = M */
        acb_calc_integrate(w, _f_stieltjes, (void *) n1, a, b, wp, tol, opt, wp);
        acb_add(v, v, w, wp);

        acb_set(a, b);
        acb_set_arb(b, M);
        arb_set(acb_imagref(b), C);  /* b = M + i C */
        acb_calc_integrate(w, _f_stieltjes, (void *) n1, a, b, wp, tol, opt, wp);
        acb_add(v, v, w, wp);

        acb_set(a, b);
        arb_set(acb_realref(b), N);  /* b = N + i C */
        acb_calc_integrate(w, _f_stieltjes, (void *) n1, a, b, wp, tol, opt, wp);
        acb_add(v, v, w, wp);

        acb_set(a, b);
        arb_zero(acb_imagref(b));    /* b = N */
        acb_calc_integrate(w, _f_stieltjes, (void *) n1, a, b, wp, tol, opt, wp);
        acb_add(v, v, w, wp);
    }

    arb_add_error_mag(acb_realref(v), bound);

    acb_const_pi(b, wp);
    acb_mul(v, v, b, wp);
    acb_div_fmpz(v, v, n1, wp);
    acb_neg(v, v);

    arb_set_round(res, acb_realref(v), prec);

    acb_clear(a);
    acb_clear(b);
    acb_clear(v);
    acb_clear(w);
    mag_clear(tol);
    mag_clear(bound);

    arb_clear(M);
    arb_clear(N);
    arb_clear(C);
}

void
acb_dirichlet_stieltjes_em(acb_t res, const fmpz_t n, const acb_t a, slong prec)
{
    if (fmpz_cmp_ui(n, 10000) > 0)
    {
        acb_indeterminate(res);
    }
    else
    {
        slong nn, wp;
        acb_ptr z;
        acb_t s;

        nn = fmpz_get_ui(n);

        acb_init(s);
        z = _acb_vec_init(nn + 1);

        wp = prec * 1.05 + 2.2*nn + 10;  /* todo: bogus */

        /* todo: we don't want to compute all the coefficients */
        acb_one(s);
        _acb_poly_zeta_cpx_series(z, s, a, 1, nn + 1, wp);

        arb_fac_ui(acb_realref(s), nn, prec + 10);
        acb_mul_arb(res, z + nn, acb_realref(s), prec);

        if (fmpz_is_odd(n))
            acb_neg(res, res);

        acb_clear(s);
        _acb_vec_clear(z, nn + 1);
    }
}

void
acb_dirichlet_stieltjes(acb_t res, const fmpz_t n, const acb_t a, slong prec)
{
    slong cutoff;

    if (acb_is_one(a) && fmpz_is_zero(n))
    {
        arb_const_euler(acb_realref(res), prec);
        arb_zero(acb_imagref(res));
        return;
    }

    if (fmpz_sgn(n) < 0)
    {
        flint_printf("stieltjes constants only defined for n >= 0");
        flint_abort();
    }

    cutoff = FLINT_MAX(100, prec / 2);
    cutoff = FLINT_MIN(cutoff, 10000);

    if (acb_is_one(a) && fmpz_cmp_ui(n, cutoff) >= 0)
    {
        fmpz_t n1;
        fmpz_init(n1);
        fmpz_add_ui(n1, n, 1);
        acb_dirichlet_stieltjes_integral(acb_realref(res), n1, prec);
        arb_zero(acb_imagref(res));
        fmpz_clear(n1);
    }
    else
    {
        acb_dirichlet_stieltjes_em(res, n, a, prec);
    }
}

