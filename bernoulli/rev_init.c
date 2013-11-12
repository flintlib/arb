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

#include "bernoulli.h"

void
bernoulli_rev_init(bernoulli_rev_t iter, ulong nmax)
{
    long j;
    fmpz_t t;
    fmprb_t x;
    int round1, round2;
    long wp;

    nmax -= (nmax % 2);
    iter->n = nmax;

    iter->alloc = 0;
    if (nmax < BERNOULLI_REV_MIN)
        return;

    iter->prec = wp = global_prec(nmax);

    iter->max_power = zeta_terms(nmax, iter->prec);
    iter->alloc = iter->max_power + 1;
    iter->powers = _fmpz_vec_init(iter->alloc);
    fmpz_init(iter->pow_error);
    fmprb_init(iter->prefactor);
    fmprb_init(iter->two_pi_squared);

    fmprb_init(x);
    fmpz_init(t);

    /* precompute powers */
    for (j = 3; j <= iter->max_power; j += 2)
    {
        fmprb_ui_pow_ui(x, j, nmax, power_prec(j, nmax, wp));
        fmprb_inv(x, x, power_prec(j, nmax, wp));
        round1 = fmpr_get_fmpz_fixed_si(t, fmprb_midref(x), -wp);
        fmpz_set(iter->powers + j, t);

        /* error: the radius, plus two roundings */
        round2 = fmpr_get_fmpz_fixed_si(t, fmprb_radref(x), -wp);
        fmpz_add_ui(t, t, (round1 != 0) + (round2 != 0));
        if (fmpz_cmp(iter->pow_error, t) < 0)
            fmpz_set(iter->pow_error, t);
    }

    /* precompute (2pi)^2 and 2*(n!)/(2pi)^n */
    fmprb_fac_ui(iter->prefactor, nmax, wp);
    fmprb_mul_2exp_si(iter->prefactor, iter->prefactor, 1);

    fmprb_const_pi(x, wp);
    fmprb_mul_2exp_si(x, x, 1);
    fmprb_mul(iter->two_pi_squared, x, x, wp);

    fmprb_pow_ui(x, iter->two_pi_squared, nmax / 2, wp);
    fmprb_div(iter->prefactor, iter->prefactor, x, wp);

    fmpz_clear(t);
    fmprb_clear(x);
}

