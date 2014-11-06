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

    Copyright (C) 2014 Fredrik Johansson

******************************************************************************/

#include "acb_modular.h"

static void
acb_mul_approx(acb_t z, acb_t tmp1, acb_t tmp2, const acb_t x, const acb_t y, long wprec, long prec)
{
    if (prec <= 1024)
    {
        acb_mul(z, x, y, wprec);
    }
    else if (x == y)
    {
        acb_set_round(tmp1, x, wprec);
        acb_mul(z, tmp1, tmp1, wprec);
    }
    else
    {
        acb_set_round(tmp1, x, wprec);
        acb_set_round(tmp2, y, wprec);
        acb_mul(z, tmp1, tmp2, wprec);
    }
}

double mag_get_log2_d_approx(const mag_t x);

#define PENTAGONAL(N) ((((N)+2)/2) * ((3*(N)+5)/2)/2)

void
acb_modular_eta_sum(acb_t eta, const acb_t q, long prec)
{
    mag_t err, qmag;
    double log2q_approx, log2term_approx;
    long e, e1, e2, k, k1, k2, N, term_prec;
    long *exponents, *aindex, *bindex;
    acb_ptr qpow;
    acb_t tmp1, tmp2;
    int q_is_real;

    acb_init(tmp1);
    acb_init(tmp2);
    mag_init(err);
    mag_init(qmag);

    q_is_real = arb_is_zero(acb_imagref(q));

    acb_get_mag(qmag, q);
    log2q_approx = mag_get_log2_d_approx(qmag);

    if (log2q_approx >= 0.0)
    {
        N = 1;
        mag_inf(err);
    }
    else  /* Pick N and compute error bound */
    {
        mag_t den;
        mag_init(den);

        N = 1;
        while (0.05 * N * N < prec)
        {
            log2term_approx = log2q_approx * PENTAGONAL(N);
            if (log2term_approx < -prec - 2)
                break;
            N++;
        }

        mag_one(den);
        mag_sub_lower(den, den, qmag);

        /* no convergence */
        if (mag_is_zero(den))
        {
            N = 1;
            mag_inf(err);
        }
        else
        {
            mag_pow_ui(err, qmag, PENTAGONAL(N));
            mag_div(err, err, den);
        }

        mag_clear(den);
    }

    exponents = flint_malloc(sizeof(long) * 3 * N);
    aindex = exponents + N;
    bindex = aindex + N;

    qpow = _acb_vec_init(N);

    acb_modular_addseq_eta(exponents, aindex, bindex, N);
    acb_set_round(qpow + 0, q, prec);

    acb_zero(eta);

    for (k = 0; k < N; k++)
    {
        e = exponents[k];

        log2term_approx = e * log2q_approx;
        term_prec = FLINT_MIN(FLINT_MAX(prec + log2term_approx + 16.0, 16.0), prec);

        if (k > 0)
        {
            k1 = aindex[k];
            k2 = bindex[k];
            e1 = exponents[k1];
            e2 = exponents[k2];

            if (e == e1 + e2)
            {
                acb_mul_approx(qpow + k, tmp1, tmp2, qpow + k1, qpow + k2, term_prec, prec);
            }
            else if (e == 2 * e1 + e2)
            {
                acb_mul_approx(qpow + k, tmp1, tmp2, qpow + k1, qpow + k1, term_prec, prec);
                acb_mul_approx(qpow + k, tmp1, tmp2, qpow + k, qpow + k2, term_prec, prec);
            }
            else
            {
                printf("exponent not in addition sequence!\n");
                abort();
            }
        }

        if (k % 4 <= 1)
            acb_sub(eta, eta, qpow + k, prec);
        else
            acb_add(eta, eta, qpow + k, prec);
    }

    acb_add_ui(eta, eta, 1, prec);

    if (q_is_real)
        arb_add_error_mag(acb_realref(eta), err);
    else
        acb_add_error_mag(eta, err);

    flint_free(exponents);
    _acb_vec_clear(qpow, N);
    acb_clear(tmp1);
    acb_clear(tmp2);
    mag_clear(err);
    mag_clear(qmag);
}

