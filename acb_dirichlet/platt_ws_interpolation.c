/*
    Copyright (C) 2019 D.H.J Polymath

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"
#include "arb_hypgeom.h"

static void
_arb_div_si_si(arb_t res, slong x, slong y, slong prec)
{
    arb_set_si(res, x);
    arb_div_si(res, res, y, prec);
}

static void
_arb_ui_pow_arb(arb_t res, ulong n, const arb_t x, slong prec)
{
    arb_t a;
    arb_init(a);
    arb_set_ui(a, n);
    arb_pow(res, a, x, prec);
    arb_clear(a);
}

/*
 * Follows the notation in https://en.wikipedia.org/wiki/Gaussian_function
 * res = a*exp(-(x-b)^2 / (2*c^2))
 * If 'a' is NULL a default coefficient of 1 is used.
 * If 'b' is NULL a default mean of 0 is used.
 */
static void
_arb_gaussian(arb_t res, const arb_t a, const arb_t b, const arb_t c,
        const arb_t x, slong prec)
{
    arb_t z;
    arb_init(z);
    if (b == NULL)
        arb_set(z, x);
    else
        arb_sub(z, x, b, prec);
    arb_div(z, z, c, prec);
    arb_sqr(z, z, prec);
    arb_mul_2exp_si(z, z, -1);
    arb_neg(z, z);
    arb_exp(z, z, prec);
    if (a != NULL)
        arb_mul(res, z, a, prec);
    arb_clear(z);
}

static void
_platt_lambda(arb_t res, const arb_t t, slong prec)
{
    acb_t pi, s, z, s1, s2;

    acb_init(pi);
    acb_init(s);
    acb_init(z);
    acb_init(s1);
    acb_init(s2);

    arb_set_d(acb_realref(s), 0.5);
    arb_set(acb_imagref(s), t);

    acb_const_pi(pi, prec);
    arb_mul_2exp_si(acb_imagref(s1), t, -1);
    acb_neg(s1, s1);
    acb_pow(s1, pi, s1, prec);

    acb_mul_2exp_si(s2, s, -1);
    acb_gamma(s2, s2, prec);

    acb_zeta(z, s, prec);
    acb_mul(z, z, s1, prec);
    acb_mul(z, z, s2, prec);
    if (!arb_contains_zero(acb_imagref(z)))
    {
        flint_abort(); 
    }
    arb_set(res, acb_realref(z));

    acb_clear(pi);
    acb_clear(s);
    acb_clear(z);
    acb_clear(s1);
    acb_clear(s2);
}

void
acb_dirichlet_platt_scaled_lambda(arb_t res, const arb_t t, slong prec)
{
    arb_t pi, lam;

    arb_init(pi);
    arb_init(lam);

    arb_const_pi(pi, prec);
    _platt_lambda(lam, t, prec);
    arb_mul(res, pi, t, prec);
    arb_mul_2exp_si(res, res, -2);
    arb_exp(res, res, prec);
    arb_mul(res, res, lam, prec);

    arb_clear(pi);
    arb_clear(lam);
}

void
acb_dirichlet_platt_scaled_lambda_vec(arb_ptr res,
        const fmpz_t T, slong A, slong B, slong prec)
{
    slong N = A*B;
    if (A < 1 || B < 1 || N % 2)
    {
        flint_printf("requires an even number of grid points\n");
        flint_abort();
    }
    else
    {
        slong i;
        arb_t t;
        arb_init(t);
        for (i = 0; i < N; i++)
        {
            slong n = i - N/2;
            _arb_div_si_si(t, n, A, prec);
            arb_add_fmpz(t, t, T, prec);
            acb_dirichlet_platt_scaled_lambda(res + i, t, prec);
        }
        arb_clear(t);
    }
}

static void
_platt_bound_C3_X(arb_t res, const arb_t t0, slong A, const arb_t H,
        slong Ns, const arb_t beta, slong prec)
{
    arb_t a, c, x;

    arb_init(a);
    arb_init(c);
    arb_init(x);

    /* res = gaussian(a=(t0+Ns/A)^beta, b=NULL, c=A*H, x=Ns) */
    _arb_div_si_si(a, Ns, A, prec);
    arb_add(a, a, t0, prec);
    arb_pow(a, a, beta, prec);
    arb_mul_si(c, H, A, prec);
    arb_set_si(x, Ns);
    _arb_gaussian(res, a, NULL, c, x, prec);

    arb_clear(a);
    arb_clear(c);
    arb_clear(x);
}

static void
_platt_bound_C3_Y(arb_t res, const arb_t t0, slong A, const arb_t H,
        slong Ns, const arb_t beta, slong prec)
{
    arb_t a, b, g, g1, g2;

    arb_init(a);
    arb_init(b);
    arb_init(g);
    arb_init(g1);
    arb_init(g2);

    /* a = 2^((2*beta - 1)/2) */
    arb_mul_2exp_si(a, beta, 1);
    arb_sub_ui(a, a, 1, prec);
    arb_mul_2exp_si(a, a, -1);
    _arb_ui_pow_arb(a, 2, a, prec);

    /* b = t0^beta */
    arb_pow(b, t0, beta, prec);

    /* g = incomplete_gamma(1/2, Ns^2 / (2*A^2*H^2)) */
    arb_set_d(g1, 0.5);
    _arb_div_si_si(g2, Ns, A, prec);
    arb_div(g2, g2, H, prec);
    arb_sqr(g2, g2, prec);
    arb_mul_2exp_si(g2, g2, -1);
    arb_hypgeom_gamma_upper(g, g1, g2, 0, prec);

    /* res = a*b*A*H*g */
    arb_mul_si(res, H, A, prec);
    arb_mul(res, res, a, prec);
    arb_mul(res, res, b, prec);
    arb_mul(res, res, g, prec);

    arb_clear(a);
    arb_clear(b);
    arb_clear(g);
    arb_clear(g1);
    arb_clear(g2);
}

static void
_platt_bound_C3_Z(arb_t res, const arb_t t0, slong A, const arb_t H,
        const arb_t beta, slong prec)
{
    arb_t a, b, g, g1, g2;

    arb_init(a);
    arb_init(b);
    arb_init(g);
    arb_init(g1);
    arb_init(g2);

    /* a = 2^((3*beta - 1)/2) */
    arb_mul_ui(a, beta, 3, prec);
    arb_sub_ui(a, a, 1, prec);
    arb_mul_2exp_si(a, a, -1);
    _arb_ui_pow_arb(a, 2, a, prec);

    /* b = H^(beta + 1) */
    arb_add_ui(b, beta, 1, prec);
    arb_pow(b, H, b, prec);

    /* g = incomplete_gamma((beta+1)/2, t0^2 / (2*H^2)) */
    arb_add_ui(g1, beta, 1, prec);
    arb_mul_2exp_si(g1, g1, -1);
    arb_div(g2, t0, H, prec);
    arb_sqr(g2, g2, prec);
    arb_mul_2exp_si(g2, g2, -1);
    arb_hypgeom_gamma_upper(g, g1, g2, 0, prec);

    /* res = a*b*A*g */
    arb_mul_si(res, g, A, prec);
    arb_mul(res, res, a, prec);
    arb_mul(res, res, b, prec);

    arb_clear(a);
    arb_clear(b);
    arb_clear(g);
    arb_clear(g1);
    arb_clear(g2);
}


/*
 * An improvement of Lemma C.3 is used, which is tighter by a factor of A.
 * https://github.com/fredrik-johansson/arb/pull/269#issuecomment-468488347
 */
void
acb_dirichlet_platt_bound_C3(arb_t res, const arb_t t0, slong A, const arb_t H,
        slong Ns, slong prec)
{
    arb_t pi, ee, beta, X, Y, Z, lhs, rhs;

    arb_init(pi);
    arb_init(ee);
    arb_init(beta);
    arb_init(X);
    arb_init(Y);
    arb_init(Z);
    arb_init(lhs);
    arb_init(rhs);

    /* t0 > exp(exp(1) is required for valid error bounds. */
    arb_const_e(ee, prec);
    arb_exp(ee, ee, prec);
    if (!arb_gt(t0, ee))
    {
        arb_zero_pm_inf(res);
        goto finish;
    }

    /* 0 < Ns <= t0*A is required for valid error bounds. */
    arb_set_si(lhs, Ns);
    arb_mul_si(rhs, t0, A, prec);
    if (!arb_is_positive(lhs) || !arb_le(lhs, rhs))
    {
        arb_zero_pm_inf(res);
        goto finish;
    }

    /* res = (X + Y + Z) * 6 / (pi * Ns) */
    arb_const_pi(pi, prec);
    acb_dirichlet_platt_beta(beta, t0, prec);
    _platt_bound_C3_X(X, t0, A, H, Ns, beta, prec);
    _platt_bound_C3_Y(Y, t0, A, H, Ns, beta, prec);
    _platt_bound_C3_Z(Z, t0, A, H, beta, prec);
    arb_add(res, X, Y, prec);
    arb_add(res, res, Z, prec);
    arb_mul_ui(res, res, 6, prec);
    arb_div(res, res, pi, prec);
    arb_div_si(res, res, Ns, prec);

finish:
    arb_clear(pi);
    arb_clear(ee);
    arb_clear(beta);
    arb_clear(X);
    arb_clear(Y);
    arb_clear(Z);
    arb_clear(lhs);
    arb_clear(rhs);
}


static void
_interpolation_helper(arb_t res, const acb_dirichlet_platt_ws_precomp_t pre,
        const arb_t t0, arb_srcptr p, const fmpz_t T, slong A, slong B,
        slong i0, slong Ns, const arb_t H, slong sigma, slong prec)
{
    arb_t dt0, dt, a, s, err, total;
    slong i;
    slong N = A*B;
    arb_init(dt0);
    arb_init(dt);
    arb_init(a);
    arb_init(s);
    arb_init(err);
    arb_init(total);
    arb_sub_fmpz(dt0, t0, T, prec + fmpz_clog_ui(T, 2));
    for (i = i0; i < i0 + 2*Ns; i++)
    {
        slong n = i - N/2;
        _arb_div_si_si(dt, n, A, prec);
        arb_sub(a, dt, dt0, prec);
        arb_mul_si(a, a, A, prec);
        arb_sinc_pi(a, a, prec);
        arb_mul(a, a, p + i, prec);
        _arb_gaussian(s, a, dt0, H, dt, prec);
        arb_add(total, total, s, prec);
    }
    acb_dirichlet_platt_bound_C3(err, t0, A, H, Ns, prec);
    arb_add_error(total, err);
    acb_dirichlet_platt_i_bound_precomp(
            err, &pre->pre_i, &pre->pre_c, t0, A, H, sigma, prec);
    arb_add_error(total, err);
    arb_set(res, total);
    arb_clear(dt0);
    arb_clear(dt);
    arb_clear(a);
    arb_clear(s);
    arb_clear(err);
    arb_clear(total);
}

void
acb_dirichlet_platt_ws_precomp_init(acb_dirichlet_platt_ws_precomp_t pre,
        slong A, const arb_t H, slong sigma, slong prec)
{
    acb_dirichlet_platt_c_precomp_init(&pre->pre_c, sigma, H, 0, prec);
    acb_dirichlet_platt_i_precomp_init(&pre->pre_i, A, H, sigma, prec);
}

void
acb_dirichlet_platt_ws_precomp_clear(acb_dirichlet_platt_ws_precomp_t pre)
{
    acb_dirichlet_platt_c_precomp_clear(&pre->pre_c);
    acb_dirichlet_platt_i_precomp_clear(&pre->pre_i);
}

void acb_dirichlet_platt_ws_interpolation_precomp(arb_t res,
    acb_dirichlet_platt_ws_precomp_t pre, const arb_t t0,
    arb_srcptr p, const fmpz_t T, slong A, slong B, slong Ns_max,
    const arb_t H, slong sigma, slong prec)
{
    slong N = A*B;
    if (A < 1 || B < 1 || N % 2)
    {
        flint_printf("requires an even number of grid points\n");
        flint_abort();
    }
    else
    {
        arb_t x, dt0, dt0A, total;
        arf_t lower_f;
        slong n, lower_n;

        arb_init(x);
        arb_init(dt0);
        arb_init(dt0A);
        arb_init(total);
        arf_init(lower_f);

        arb_sub_fmpz(dt0, t0, T, prec + fmpz_clog_ui(T, 2));
        arb_mul_si(dt0A, dt0, A, prec);
        arb_get_lbound_arf(lower_f, dt0A, prec);
        lower_n = arf_get_si(lower_f, ARF_RND_FLOOR);

        /*
         * More than one iteration is needed only when the set of
         * supporting points for interpolation is uncertain.
         */
        for (n = lower_n; n == lower_n || arb_contains_si(dt0A, n); n++)
        {
            slong nlow = N/2 + n + 1;
            slong nhigh = N/2 - n - 1;
            slong Ns = FLINT_MIN(Ns_max, FLINT_MIN(nlow, nhigh));
            if (Ns < 1)
            {
                arb_zero_pm_inf(total);
            }
            else
            {
                slong i0 = N/2 + n - (Ns - 1);
                _interpolation_helper(
                        x, pre, t0, p, T, A, B, i0, Ns, H, sigma, prec);
                if (n == lower_n)
                {
                    arb_set(total, x);
                }
                else
                {
                    arb_union(total, total, x, prec);
                }
            }
        }

        arb_set(res, total);

        arb_clear(x);
        arb_clear(dt0);
        arb_clear(dt0A);
        arb_clear(total);
        arf_clear(lower_f);
    }
}

void
acb_dirichlet_platt_ws_interpolation(arb_t res, const arb_t t0,
        arb_srcptr p, const fmpz_t T, slong A, slong B,
        slong Ns_max, const arb_t H, slong sigma, slong prec)
{
    acb_dirichlet_platt_ws_precomp_t pre;
    acb_dirichlet_platt_ws_precomp_init(pre, A, H, sigma, prec);
    acb_dirichlet_platt_ws_interpolation_precomp(
            res, pre, t0, p, T, A, B, Ns_max, H, sigma, prec);
    acb_dirichlet_platt_ws_precomp_clear(pre);
}
