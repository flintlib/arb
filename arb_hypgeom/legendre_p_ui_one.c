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

static void PP(fmpz_t r, ulong n, ulong k, ulong prime)
{
    if (n < BIG && k < BIG)
    {
        fmpz_set_ui(r, (n - k + 1 - prime) * (n + k + prime));
    }
    else
    {
        fmpz_set_ui(r, n - k + 1 - prime);
        fmpz_mul_ui(r, r, n + k + prime);
    }
}

static void QQ(fmpz_t r, ulong n, ulong k, ulong prime)
{
    if (k < BIG)
    {
        fmpz_set_ui(r, k * (k + prime));
    }
    else
    {
        fmpz_set_ui(r, k);
        fmpz_mul_ui(r, r, k + prime);
    }
}

static void
sum_rs_inner(arb_t s, arb_srcptr xpow, slong m, ulong n, slong K, ulong prime, slong prec)
{
    slong j, k, khi, klo, u, r;
    fmpz_t cc, dd;

    fmpz_init(cc);
    fmpz_init(dd);

    arb_zero(s);
    fmpz_zero(cc);

    k = K - 1;
    while (k >= 1)
    {
        u = FLINT_MIN(5, k);
        khi = k;
        klo = k - u + 1;

        PP(cc, n, khi, prime);
        for (j = klo; j < khi; j++)
        {
            PP(dd, n, j, prime);
            fmpz_mul(cc, cc, dd);
        }

        while (k >= klo)
        {
            r = k % m;
            if (k == khi)
            {
                arb_add(s, s, xpow + r, prec);
                arb_mul_fmpz(s, s, cc, prec);
            }
            else if (r == 0)
                arb_add_fmpz(s, s, cc, prec);
            else
                arb_addmul_fmpz(s, xpow + r, cc, prec);

            if (r == 0 && k != 0)
                arb_mul(s, s, xpow + m, prec);

            PP(dd, n, k, prime);
            fmpz_divexact(cc, cc, dd);
            QQ(dd, n, k, prime);
            fmpz_mul(cc, cc, dd);
            k--;
        }

        arb_div_fmpz(s, s, cc, prec);
    }

    fmpz_clear(cc);
    fmpz_clear(dd);
}

void
arb_hypgeom_legendre_p_ui_one(arb_t res, arb_t res_prime, ulong n,
                        const arb_t x, slong K, slong prec)
{
    arb_t s, v;
    arb_ptr xpow;
    slong m;
    mag_t u, a, t, err;

    if (n == 0)
    {
        if (res != NULL) arb_one(res);
        if (res_prime != NULL) arb_zero(res_prime);
        return;
    }

    /* overflow protection */
    if (n > UWORD_MAX / 4)
    {
        if (res != NULL) arb_indeterminate(res);
        if (res_prime != NULL) arb_indeterminate(res_prime);
    }

    arb_init(v);
    arb_init(s);
    mag_init(u);
    mag_init(a);
    mag_init(t);
    mag_init(err);

    K = FLINT_MIN(K, n + 1);

    if (res != NULL && res_prime != NULL)
        m = n_sqrt(2 * K);
    else
        m = n_sqrt(K);

    xpow = _arb_vec_init(m + 1);

    arb_sub_ui(v, x, 1, prec);
    arb_mul_2exp_si(v, v, -1);
    _arb_vec_set_powers(xpow, v, m + 1, prec);

    /* truncating */
    if (K < n + 1)
    {
        arb_get_mag(u, v);
        mag_mul_ui(t, u, n - K);
        mag_mul_ui(t, t, n + K + 1);
        mag_div_ui(t, t, K + 1);
        mag_div_ui(t, t, K + 1);
        mag_geom_series(t, t, 0);
        mag_pow_ui(u, u, K);
        mag_mul(u, u, t);
    }

    if (res != NULL)
    {
        sum_rs_inner(s, xpow, m, n, K, 0, prec);
        arb_add_ui(res, s, 1, prec);

        if (K < n + 1)
        {
            mag_set(err, u);
            mag_bin_uiui(t, n, K);
            mag_mul(err, err, t);
            mag_bin_uiui(t, n + K, K);
            mag_mul(err, err, t);
            arb_add_error_mag(res, err);
        }
    }

    if (res_prime != NULL)
    {
        K = FLINT_MIN(K, n);
        sum_rs_inner(s, xpow, m, n, K, 1, prec);
        arb_add_ui(res_prime, s, 1, prec);
        arb_mul_ui(res_prime, res_prime, n, prec);
        arb_mul_ui(res_prime, res_prime, n + 1, prec);
        arb_mul_2exp_si(res_prime, res_prime, -1);

        /* truncating */
        if (K < n)
        {
            mag_set(err, u);
            mag_bin_uiui(t, n, K + 1);
            mag_mul(err, err, t);
            mag_bin_uiui(t, n + K + 1, K + 1);
            mag_mul(err, err, t);
            mag_mul_ui(err, err, n);
            arb_add_error_mag(res_prime, err);
        }

    }

    _arb_vec_clear(xpow, m + 1);
    arb_clear(s);
    arb_clear(v);
    mag_clear(u);
    mag_clear(a);
    mag_clear(t);
    mag_clear(err);
}

