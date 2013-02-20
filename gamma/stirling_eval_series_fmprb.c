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

void
gamma_stirling_bound_remainder(fmpr_t err, const fmprb_t z, long n)
{
    if (fmprb_contains_nonpositive(z))
    {
        fmpr_pos_inf(err);
    }
    else
    {
        fmpr_t t;
        fmprb_t b;

        fmpr_init(t);
        fmprb_init(b);

        fmpr_sub(t, fmprb_midref(z), fmprb_radref(z), FMPRB_RAD_PREC, FMPR_RND_FLOOR);

        if (fmpr_sgn(t) <= 0)
        {
            fmpr_pos_inf(err);
        }
        else
        {
            fmpr_ui_div(t, 1, t, FMPRB_RAD_PREC, FMPR_RND_UP);
            fmpr_pow_sloppy_ui(t, t, 2 * n - 1, FMPRB_RAD_PREC, FMPR_RND_UP);

            gamma_stirling_coeff(b, n, FMPRB_RAD_PREC);

            fmprb_get_abs_ubound_fmpr(err, b, FMPRB_RAD_PREC);
            fmpr_mul(err, err, t, FMPRB_RAD_PREC, FMPR_RND_UP);
        }

        fmpr_clear(t);
        fmprb_clear(b);
    }
}

void
gamma_stirling_eval_series_fmprb(fmprb_t s, const fmprb_t z, long nterms, long prec)
{
    fmprb_t t, u, b, w, v;
    fmpr_t err;

    long k, term_prec;
    double z_mag, term_mag;

    fmprb_init(t);
    fmprb_init(u);
    fmprb_init(b);
    fmprb_init(w);
    fmprb_init(v);

    fmprb_log(w, z, prec);

    bernoulli_cache_compute(2 * (nterms + 1));

    nterms = FLINT_MAX(nterms, 1);

    fmprb_zero(s);
    if (nterms > 1)
    {
        fmprb_ui_div(t, 1UL, z, prec);
        fmprb_mul(u, t, t, prec);
        z_mag = fmpr_get_d(fmprb_midref(w), FMPR_RND_UP) * 1.44269504088896;

        for (k = nterms - 1; k >= 1; k--)
        {
            term_mag = bernoulli_bound_2exp_si(2 * k);
            term_mag -= (2 * k - 1) * z_mag;
            term_prec = prec + term_mag;
            term_prec = FLINT_MIN(term_prec, prec);
            term_prec = FLINT_MAX(term_prec, 10);

            gamma_stirling_coeff(b, k, term_prec);

            if (prec > 2000)
            {
                fmprb_set_round(v, u, term_prec);
                fmprb_mul(s, s, v, term_prec);
            }
            else
            {
                fmprb_mul(s, s, u, term_prec);
            }

            fmprb_add(s, s, b, term_prec);
        }

        fmprb_mul(s, s, t, prec);
    }

    /* remainder bound */
    fmpr_init(err);
    gamma_stirling_bound_remainder(err, z, nterms);
    fmprb_add_error_fmpr(s, err);
    fmpr_clear(err);

    /* (z-0.5)*log(z) - z + log(2*pi)/2 */
    fmprb_set_ui(t, 1);
    fmprb_mul_2exp_si(t, t, -1);
    fmprb_sub(t, z, t, prec);
    fmprb_mul(t, w, t, prec);

    fmprb_add(s, s, t, prec);
    fmprb_sub(s, s, z, prec);

    fmprb_const_log_sqrt2pi(t, prec);
    fmprb_add(s, s, t, prec);

    fmprb_clear(t);
    fmprb_clear(u);
    fmprb_clear(b);
    fmprb_clear(w);
    fmprb_clear(v);
}

