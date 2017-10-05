/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_hypgeom.h"

#define BIG (1 << (FLINT_BITS / 2 - 2))

static void PP(fmpz_t res, ulong n, slong k)
{
    /* (2 * k - 1) * (2 * k - 1) */
    if (k < BIG)
    {
        fmpz_set_ui(res, (2 * k - 1) * (2 * k - 1));
    }
    else
    {
        fmpz_set_ui(res, 2 * k - 1);
        fmpz_mul_ui(res, res, 2 * k - 1);
    }
}

static void QQ(fmpz_t res, ulong n, slong k)
{
    if (k < BIG && n < BIG)
    {
        fmpz_set_ui(res, k * (2 * k + 2 * n + 1));
    }
    else
    {
        fmpz_set_ui(res, n);
        fmpz_mul_2exp(res, res, 1);
        fmpz_add_ui(res, res, 2 * k + 1);
        fmpz_mul_ui(res, res, k);
    }
}

static void
asymp_series(acb_t res, ulong n, const acb_t x, slong K, slong prec)
{
    slong j, k, khi, klo, m, u, r;
    acb_t s;
    acb_ptr xpow;
    fmpz_t cc, dd;

    m = n_sqrt(K);

    acb_init(s);
    fmpz_init(cc);
    fmpz_init(dd);
    xpow = _acb_vec_init(m + 1);

    _acb_vec_set_powers(xpow, x, m + 1, prec);

    k = K - 1;
    while (k >= 1)
    {
        u = FLINT_MIN(5, k);
        khi = k;
        klo = k - u + 1;

        PP(cc, n, khi);
        for (j = klo; j < khi; j++)
        {
            PP(dd, n, j);
            fmpz_mul(cc, cc, dd);
        }

        while (k >= klo)
        {
            r = k % m;

            if (k == khi)
            {
                acb_add(s, s, xpow + r, prec);
                acb_mul_fmpz(s, s, cc, prec);
            }
            else if (r == 0)
            {
                acb_add_fmpz(s, s, cc, prec);
            }
            else
            {
                acb_addmul_fmpz(s, xpow + r, cc, prec);
            }

            if (r == 0 && k != 0)
                acb_mul(s, s, xpow + m, prec);

            PP(dd, n, k);
            fmpz_divexact(cc, cc, dd);
            QQ(dd, n, k);
            fmpz_mul(cc, cc, dd);

            k--;
        }

        acb_div_fmpz(s, s, cc, prec);
    }

    acb_add_ui(res, s, 1, prec);

    _acb_vec_clear(xpow, m + 1);
    acb_clear(s);
    fmpz_clear(cc);
    fmpz_clear(dd);
}


/* error: 0.50795 / (y^K sqrt(ny)) * [K! n! / (2^K (n+K)!)] */
/*                                   [K! / (2n)^K]          */
static void
_arb_hypgeom_legendre_p_ui_asymp_error(mag_t res, ulong n, const mag_t y, slong K)
{
    mag_t t, u;
    mag_init(t);
    mag_init(u);

    /* t = K! / (y^K sqrt(ny)) */
    mag_mul_ui_lower(t, y, n);
    mag_sqrt_lower(t, t);
    mag_pow_ui_lower(u, y, K);
    mag_mul_lower(t, t, u);
    mag_fac_ui(u, K);
    mag_div(t, u, t);

    if (K < n / 16)
    {
        /* (2n)^K */
        mag_set_ui_lower(u, n);
        mag_mul_2exp_si(u, u, 1);
        mag_pow_ui_lower(u, u, K);
        mag_div(t, t, u);
    }
    else
    {
        /* n! */
        mag_fac_ui(u, n);
        mag_mul(t, t, u);
        /* (n+K)! */
        mag_rfac_ui(u, n + K);
        mag_mul(t, t, u);
        /* 2^K */
        mag_mul_2exp_si(t, t, -K);
    }

    mag_mul_ui(t, t, 131);
    mag_mul_2exp_si(t, t, -8);

    mag_set(res, t);

    mag_clear(t);
    mag_clear(u);
}

int
arb_abs_le_ui(const arb_t x, ulong n)
{
    arf_struct u[3];
    arf_t t;
    int res;

    if (!arb_is_finite(x) || arf_cmpabs_ui(arb_midref(x), n) > 0)
        return 0;

    if (arb_is_exact(x))
        return 1;

    if (arf_sgn(arb_midref(x)) >= 0)
        arf_init_set_shallow(u + 0, arb_midref(x));
    else
        arf_init_neg_shallow(u + 0, arb_midref(x));

    arf_init_set_mag_shallow(u + 1, arb_radref(x));

    arf_init(u + 2);   /* no need to free */
    arf_set_ui(u + 2, n);
    arf_neg(u + 2, u + 2);

    arf_init(t);
    arf_sum(t, u, 3, MAG_BITS, ARF_RND_DOWN);
    res = (arf_sgn(t) < 0);
    arf_clear(t);

    return res;
}

void
arb_hypgeom_legendre_p_ui_asymp(arb_t res, ulong n, const arb_t x, slong K, slong prec)
{
    arb_t y, t, u;
    acb_t w, s, z;
    fmpz_t e;
    mag_t err;

    if (!arb_abs_le_ui(x, 1))
    {
        arb_indeterminate(res);
        return;
    }

    K = FLINT_MAX(K, 1);

    arb_init(y);
    arb_init(t);
    arb_init(u);
    acb_init(w);
    acb_init(s);
    acb_init(z);
    mag_init(err);
    fmpz_init(e);

    /* u = n + 1/2 */
    arb_set_d(u, 0.5);
    arb_add_ui(u, u, n, prec);

    /* y = sqrt(1-x^2) */
    arb_one(y);
    arb_submul(y, x, x, 2 * prec);
    arb_sqrt(y, y, prec);

    arb_get_mag_lower(err, y);
    _arb_hypgeom_legendre_p_ui_asymp_error(err, n, err, K);

    /* w = 1 - (x/y)i */
    arb_one(acb_realref(w));
    arb_div(acb_imagref(w), x, y, prec);
    arb_neg(acb_imagref(w), acb_imagref(w));

    /* z = (x + yi)^(n+0.5) * (1-i) */
    arb_set(acb_realref(z), x);
    arb_set(acb_imagref(z), y);
    acb_pow_arb(z, z, u, prec);
    arb_one(acb_realref(s));
    arb_set_si(acb_imagref(s), -1);
    acb_mul(z, z, s, prec);

    /* main series */
    acb_mul_2exp_si(s, w, -2);
    asymp_series(s, n, s, K, prec);

    /* we will use Re(z * s) */
    acb_mul(z, z, s, prec);

    /* prefactor: t = 4^n / (pi * sqrt(y) * (n+0.5) binomial(2n,n)) */
    arb_hypgeom_central_bin_ui(t, n, prec);

    arb_mul(t, t, u, prec);
    arb_sqrt(u, y, prec);
    arb_mul(t, t, u, prec);
    arb_const_pi(u, prec);
    arb_mul(t, t, u, prec);

    arb_div(res, acb_realref(z), t, prec);

    fmpz_set_ui(e, n);
    arb_mul_2exp_fmpz(res, res, e);
    arb_mul_2exp_fmpz(res, res, e);

    arb_add_error_mag(res, err);

    arb_clear(y);
    arb_clear(t);
    arb_clear(u);
    acb_clear(w);
    acb_clear(s);
    acb_clear(z);
    mag_clear(err);
    fmpz_clear(e);
}

