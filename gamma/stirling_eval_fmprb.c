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
gamma_stirling_eval_fmprb(fmprb_t s, const fmprb_t z, long nterms, int digamma, long prec)
{
    fmprb_t b, t, logz, zinv, zinv2;
    fmpr_t err;

    long k, term_prec;
    double z_mag, term_mag;

    fmprb_init(b);
    fmprb_init(t);
    fmprb_init(logz);
    fmprb_init(zinv);
    fmprb_init(zinv2);

    fmprb_log(logz, z, prec);
    fmprb_ui_div(zinv, 1UL, z, prec);

    nterms = FLINT_MAX(nterms, 1);

    fmprb_zero(s);

    if (nterms > 1)
    {
        fmprb_mul(zinv2, zinv, zinv, prec);

        z_mag = fmpr_get_d(fmprb_midref(logz), FMPR_RND_UP) * 1.44269504088896;

        for (k = nterms - 1; k >= 1; k--)
        {
            term_mag = bernoulli_bound_2exp_si(2 * k);
            term_mag -= (2 * k - 1) * z_mag;
            term_prec = prec + term_mag;
            term_prec = FLINT_MIN(term_prec, prec);
            term_prec = FLINT_MAX(term_prec, 10);

            if (prec > 2000)
            {
                fmprb_set_round(t, zinv2, term_prec);
                fmprb_mul(s, s, t, term_prec);
            }
            else
                fmprb_mul(s, s, zinv2, term_prec);

            gamma_stirling_coeff(b, k, digamma, term_prec);
            fmprb_add(s, s, b, term_prec);
        }

        if (digamma)
            fmprb_mul(s, s, zinv2, prec);
        else
            fmprb_mul(s, s, zinv, prec);
    }

    /* remainder bound */
    fmpr_init(err);
    gamma_stirling_bound_fmprb(err, z, digamma ? 1 : 0, 1, nterms);
    fmprb_add_error_fmpr(s, err);
    fmpr_clear(err);

    if (digamma)
    {
        fmprb_neg(s, s);
        fmprb_mul_2exp_si(zinv, zinv, -1);
        fmprb_sub(s, s, zinv, prec);
        fmprb_add(s, s, logz, prec);
    }
    else
    {
        /* (z-0.5)*log(z) - z + log(2*pi)/2 */
        fmprb_one(t);
        fmprb_mul_2exp_si(t, t, -1);
        fmprb_sub(t, z, t, prec);
        fmprb_mul(t, logz, t, prec);
        fmprb_add(s, s, t, prec);
        fmprb_sub(s, s, z, prec);
        fmprb_const_log_sqrt2pi(t, prec);
        fmprb_add(s, s, t, prec);
    }

    fmprb_clear(t);
    fmprb_clear(b);
    fmprb_clear(zinv);
    fmprb_clear(zinv2);
    fmprb_clear(logz);
}

