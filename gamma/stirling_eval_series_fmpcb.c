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
#include "gamma.h"
#include "bernoulli.h"

/*
B_(2n) / (2n (2n-1)) / |z|^(2n-1) * (1/cos(0.5*arg(z))^(2n))
*/
void
gamma_stirling_bound_remainder_fmpcb(fmpr_t err, const fmpcb_t z, int digamma, long n)
{
    fmpr_t t;
    fmprb_t b;

    if (fmprb_contains_zero(fmpcb_imagref(z)) &&
        fmprb_contains_nonpositive(fmpcb_realref(z)))
    {
        fmpr_pos_inf(err);
        return;
    }

    /* bound 1 / |z|^(2n-1) */
    fmpcb_get_abs_lbound_fmpr(err, z, FMPRB_RAD_PREC);

    if (fmpr_is_zero(err))
    {
        fmpr_pos_inf(err);
        return;
    }
    fmpr_ui_div(err, 1, err, FMPRB_RAD_PREC, FMPR_RND_UP);

    if (digamma)
        fmpr_pow_sloppy_ui(err, err, 2 * n, FMPRB_RAD_PREC, FMPR_RND_UP);
    else
        fmpr_pow_sloppy_ui(err, err, 2 * n - 1, FMPRB_RAD_PREC, FMPR_RND_UP);

    /* bound coefficient */
    fmprb_init(b);
    fmpr_init(t);

    gamma_stirling_coeff(b, n, digamma, FMPRB_RAD_PREC);
    fmprb_get_abs_ubound_fmpr(t, b, FMPRB_RAD_PREC);
    fmpr_mul(err, err, t, FMPRB_RAD_PREC, FMPR_RND_UP);

    /* bound 1/cos(0.5*arg(z))^(2n) */
    fmpcb_arg(b, z, FMPRB_RAD_PREC);
    fmprb_mul_2exp_si(b, b, -1);
    fmprb_cos(b, b, FMPRB_RAD_PREC);
    fmprb_get_abs_lbound_fmpr(t, b, FMPRB_RAD_PREC);
    fmpr_ui_div(t, 1, t, FMPRB_RAD_PREC, FMPR_RND_UP);

    if (digamma)
        fmpr_pow_sloppy_ui(t, t, 2 * n + 1, FMPRB_RAD_PREC, FMPR_RND_UP);
    else
        fmpr_pow_sloppy_ui(t, t, 2 * n, FMPRB_RAD_PREC, FMPR_RND_UP);

    fmpr_mul(err, err, t, FMPRB_RAD_PREC, FMPR_RND_UP);

    fmprb_clear(b);
    fmpr_clear(t);
}

void
gamma_stirling_eval_series_fmpcb(fmpcb_t s, const fmpcb_t z, long nterms, int digamma, long prec)
{
    fmpcb_t t, logz, zinv, zinv2;
    fmprb_t b;
    fmpr_t err;

    long k, term_prec;
    double z_mag, term_mag;

    fmpcb_init(t);
    fmpcb_init(logz);
    fmpcb_init(zinv);
    fmpcb_init(zinv2);
    fmprb_init(b);

    fmpcb_log(logz, z, prec);
    fmpcb_inv(zinv, z, prec);

    nterms = FLINT_MAX(nterms, 1);

    fmpcb_zero(s);
    if (nterms > 1)
    {
        fmpcb_mul(zinv2, zinv, zinv, prec);

        z_mag = fmpr_get_d(fmprb_midref(fmpcb_realref(logz)), FMPR_RND_UP) * 1.44269504088896;

        for (k = nterms - 1; k >= 1; k--)
        {
            term_mag = bernoulli_bound_2exp_si(2 * k);
            term_mag -= (2 * k - 1) * z_mag;
            term_prec = prec + term_mag;
            term_prec = FLINT_MIN(term_prec, prec);
            term_prec = FLINT_MAX(term_prec, 10);

            gamma_stirling_coeff(b, k, digamma, term_prec);

            if (prec > 2000)
            {
                fmpcb_set_round(t, zinv2, term_prec);
                fmpcb_mul(s, s, t, term_prec);
            }
            else
                fmpcb_mul(s, s, zinv2, term_prec);

            fmprb_add(fmpcb_realref(s), fmpcb_realref(s), b, term_prec);
        }

        if (digamma)
            fmpcb_mul(s, s, zinv2, prec);
        else
            fmpcb_mul(s, s, zinv, prec);
    }

    /* remainder bound */
    fmpr_init(err);
    gamma_stirling_bound_remainder_fmpcb(err, z, digamma, nterms);
    fmprb_add_error_fmpr(fmpcb_realref(s), err);
    fmprb_add_error_fmpr(fmpcb_imagref(s), err);
    fmpr_clear(err);

    if (digamma)
    {
        fmpcb_neg(s, s);
        fmpcb_mul_2exp_si(zinv, zinv, -1);
        fmpcb_sub(s, s, zinv, prec);
        fmpcb_add(s, s, logz, prec);
    }
    else
    {
        /* (z-0.5)*log(z) - z + log(2*pi)/2 */
        fmprb_one(b);
        fmprb_mul_2exp_si(b, b, -1);
        fmprb_set(fmpcb_imagref(t), fmpcb_imagref(z));
        fmprb_sub(fmpcb_realref(t), fmpcb_realref(z), b, prec);
        fmpcb_mul(t, logz, t, prec);
        fmpcb_add(s, s, t, prec);
        fmpcb_sub(s, s, z, prec);
        fmprb_const_log_sqrt2pi(b, prec);
        fmprb_add(fmpcb_realref(s), fmpcb_realref(s), b, prec);
    }

    fmpcb_clear(t);
    fmpcb_clear(logz);
    fmpcb_clear(zinv);
    fmpcb_clear(zinv2);
    fmprb_clear(b);
}

