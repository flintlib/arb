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

void mag_sub_lower(mag_t z, const mag_t x, const mag_t y);

double
mag_get_log2_d_approx(const mag_t x)
{
    if (mag_is_zero(x))
    {
        return COEFF_MIN;
    }
    else if (mag_is_inf(x))
    {
        return COEFF_MAX;
    }
    else if (COEFF_IS_MPZ(MAG_EXP(x)))
    {
        if (fmpz_sgn(MAG_EXPREF(x)) < 0)
            return COEFF_MIN;
        else
            return COEFF_MAX;
    }
    else
    {
        long e = MAG_EXP(x);

        if (e < -20 || e > 20)
            return e;
        else
            return e + 1.4426950408889634074 *
                mag_d_log_upper_bound(MAG_MAN(x) * ldexp(1.0, -MAG_BITS));
    }
}

void
acb_modular_theta_1234_sum(acb_t theta1, acb_t theta2,
        acb_t theta3, acb_t theta4,
    const acb_t w, int w_is_unit, const acb_t q, long prec)
{
    mag_t err, qmag, wmag, vmag;
    double log2q_approx, log2w_approx, log2term_approx;
    long e, e1, e2, k, k1, k2, N, WN, term_prec;
    long *exponents, *aindex, *bindex;
    acb_ptr qpow, wpow, vpow;
    acb_t tmp1, tmp2, v;
    int q_is_real, w_is_one;

    q_is_real = arb_is_zero(acb_imagref(q));
    w_is_one = acb_is_one(w);

    acb_init(tmp1);
    acb_init(tmp2);
    acb_init(v);
    mag_init(err);

    mag_init(qmag);
    mag_init(wmag);
    mag_init(vmag);

    if (w_is_one)
        acb_one(v);
    else if (w_is_unit)
        acb_conj(v, w);
    else
        acb_inv(v, w, prec);

    acb_get_mag(qmag, q);
    log2q_approx = mag_get_log2_d_approx(qmag);

    if (w_is_unit)
    {
        mag_one(wmag);
        mag_one(vmag);
        log2w_approx = 0.0;
    }
    else
    {
        acb_get_mag(wmag, w);
        acb_get_mag(vmag, v);
        mag_max(wmag, wmag, vmag);
        log2w_approx = mag_get_log2_d_approx(wmag);
    }

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
            log2term_approx = log2q_approx * ((N+2)*(N+2)/4) + (N+2)*log2w_approx;
            if (log2term_approx < -prec - 2)
                break;
            N++;
        }

        if (w_is_unit)
        {
            mag_one(den);
            mag_sub_lower(den, den, qmag);  /* 1 - |q| is good enough */
        }
        else  /* denominator: 1 - |q|^(floor((N+1)/2)+1) * max(|w|,1/|w|) */
        {
            mag_pow_ui(err, qmag, (N + 1) / 2 + 1);
            mag_mul(err, err, wmag);
            mag_one(den);
            mag_sub_lower(den, den, err);
        }

        /* no convergence */
        if (mag_is_zero(den))
        {
            N = 1;
            mag_inf(err);
        }
        else if (w_is_unit)
        {
            mag_pow_ui(err, qmag, ((N + 2) * (N + 2)) / 4);
            mag_div(err, err, den);
            mag_mul_2exp_si(err, err, 1);
        }
        else
        {
            mag_pow_ui(err, qmag, ((N + 2) * (N + 2)) / 4);
            mag_pow_ui(vmag, wmag, N + 2);
            mag_mul(err, err, vmag);
            mag_div(err, err, den);
            mag_mul_2exp_si(err, err, 1);
        }

        mag_clear(den);
    }

    exponents = flint_malloc(sizeof(long) * 3 * N);
    aindex = exponents + N;
    bindex = aindex + N;

    qpow = _acb_vec_init(N);

    acb_modular_addseq_theta(exponents, aindex, bindex, N);
    acb_set_round(qpow + 0, q, prec);

    acb_zero(theta1);
    acb_zero(theta2);
    acb_zero(theta3);
    acb_zero(theta4);

    WN = (N + 3) / 2;

    /* compute powers of w^2 and v = 1/w^2 */
    if (!w_is_one)
    {
        wpow = _acb_vec_init(WN);
        vpow = _acb_vec_init(WN + 1);

        acb_mul(tmp1, w, w, prec);
        acb_mul(tmp2, v, v, prec);

        _acb_vec_set_powers(wpow, tmp1, WN, prec);
        _acb_vec_set_powers(vpow, tmp2, WN + 1, prec);
    }
    else
    {
        wpow = vpow = NULL;
    }

    for (k = 0; k < N; k++)
    {
        e = exponents[k];

        log2term_approx = e * log2q_approx + (k+2) * log2w_approx;
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

        if (w_is_one)
        {
            if (k % 2 == 0)
            {
                acb_add(theta3, theta3, qpow + k, prec);

                if (k % 4 == 0)
                    acb_sub(theta4, theta4, qpow + k, prec);
                else
                    acb_add(theta4, theta4, qpow + k, prec);
            }
            else
            {
                acb_add(theta2, theta2, qpow + k, prec);
            }
        }
        else
        {
            if (k % 2 == 0)
            {
                acb_add(tmp1, wpow + k / 2 + 1, vpow + k / 2 + 1, term_prec);
                acb_mul(tmp1, qpow + k, tmp1, term_prec);

                acb_add(theta3, theta3, tmp1, prec);

                if (k % 4 == 0)
                    acb_sub(theta4, theta4, tmp1, prec);
                else
                    acb_add(theta4, theta4, tmp1, prec);
            }
            else
            {
                if (k / 2 + 1 > WN - 1)
                    abort();
                if (k / 2 + 2 > WN + 1 - 1)
                    abort();

                acb_add(tmp1, wpow + k / 2 + 1, vpow + k / 2 + 2, term_prec);
                acb_mul(tmp1, qpow + k, tmp1, term_prec);
                acb_add(theta2, theta2, tmp1, prec);

                acb_sub(tmp1, wpow + k / 2 + 1, vpow + k / 2 + 2, term_prec);
                acb_mul(tmp1, qpow + k, tmp1, term_prec);
                if (k % 4 == 1)
                    acb_sub(theta1, theta1, tmp1, prec);
                else
                    acb_add(theta1, theta1, tmp1, prec);
            }
        }
    }

    if (w_is_one)
    {
        acb_mul_2exp_si(theta2, theta2, 1);
        acb_mul_2exp_si(theta3, theta3, 1);
        acb_mul_2exp_si(theta4, theta4, 1);

        acb_add_ui(theta2, theta2, 2, prec);
        acb_add_ui(theta3, theta3, 1, prec);
        acb_add_ui(theta4, theta4, 1, prec);
    }
    else
    {
        /* w * [(1 - w^-2) + series] */
        acb_sub(theta1, theta1, vpow + 1, prec);
        acb_mul(theta1, theta1, w, prec);
        acb_add(theta1, theta1, w, prec);

        /* multiply by -i */
        acb_mul_onei(theta1, theta1);
        acb_neg(theta1, theta1);

        /* w * [(1 + w^-2) + series] */
        acb_add(theta2, theta2, vpow + 1, prec);
        acb_mul(theta2, theta2, w, prec);
        acb_add(theta2, theta2, w, prec);

        acb_add_ui(theta3, theta3, 1, prec);
        acb_add_ui(theta4, theta4, 1, prec);
    }

    if (q_is_real && w_is_unit)    /* result must be real */
    {
        arb_add_error_mag(acb_realref(theta1), err);
        arb_add_error_mag(acb_realref(theta2), err);
        arb_add_error_mag(acb_realref(theta3), err);
        arb_add_error_mag(acb_realref(theta4), err);
    }
    else
    {
        acb_add_error_mag(theta1, err);
        acb_add_error_mag(theta2, err);
        acb_add_error_mag(theta3, err);
        acb_add_error_mag(theta4, err);
    }

    if (!w_is_one)
    {
        _acb_vec_clear(wpow, WN);
        _acb_vec_clear(vpow, WN + 1);
    }

    flint_free(exponents);
    _acb_vec_clear(qpow, N);
    acb_clear(tmp1);
    acb_clear(tmp2);
    acb_clear(v);
    mag_clear(err);
    mag_clear(qmag);
    mag_clear(wmag);
    mag_clear(vmag);
}

