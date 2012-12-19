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
bernoulli_rev_next(fmpz_t numer, fmpz_t denom, bernoulli_rev_t iter)
{
    ulong n;
    long j, wp;
    fmpz_t sum;
    fmpr_t err;
    fmprb_t z, h;

    n = iter->n;
    wp = iter->prec;

    if (n < bernoulli_rev_MIN)
    {
        _arith_bernoulli_number(numer, denom, n);
        if (n != 0)
            iter->n -= 2;
        return;
    }

    fmpz_init(sum);
    fmpr_init(err);
    fmprb_init(z);
    fmprb_init(h);

    /* add all odd powers */
    fmpz_zero(sum);
    for (j = iter->max_power; j >= 3; j -= 2)
        fmpz_add(sum, sum, iter->powers + j);
    fmprb_set_fmpz(z, sum);

    /* bound numerical error from the powers */
    fmpz_mul_ui(sum, iter->pow_error, iter->max_power / 2);
    fmpr_set_fmpz(err, sum);
    fmprb_add_error_fmpr(z, err);
    fmprb_mul_2exp_si(z, z, -wp);

    fmprb_add_ui(z, z, 1, wp);

    /* add truncation error: sum_{k > N} 1/k^n <= 1/N^(i-1) */
    fmpr_set_ui(err, iter->max_power);
    fmpr_pow_sloppy_ui(err, err, n - 1, FMPRB_RAD_PREC, FMPR_RND_DOWN);
    fmpr_ui_div(err, 1, err, FMPRB_RAD_PREC, FMPR_RND_UP);
    fmprb_add_error_fmpr(z, err);

    /* convert zeta to Bernoulli number */
    fmprb_div_2expm1_ui(h, z, n, wp);
    fmprb_add(z, z, h, wp);
    fmprb_mul(z, z, iter->prefactor, wp);
    arith_bernoulli_number_denom(denom, n);
    fmprb_mul_fmpz(z, z, denom, wp);
    if (n % 4 == 0)
        fmprb_neg(z, z);

    /* printf("%ld: ", n); fmprb_printd(z, 5); printf("\n"); */

    if (!fmprb_get_unique_fmpz(numer, z))
    {
        printf("warning: insufficient precision for B_%ld\n", n);
        _arith_bernoulli_number(numer, denom, n);
    }

    /* update prefactor */
    if (n > 0)
    {
        fmprb_mul(iter->prefactor, iter->prefactor, iter->two_pi_squared, wp);
        fmprb_div_ui(iter->prefactor, iter->prefactor, n, wp);
        fmprb_div_ui(iter->prefactor, iter->prefactor, n - 1, wp);
    }

    /* update powers */
    for (j = 3; j <= iter->max_power; j += 2)
        fmpz_mul2_uiui(iter->powers + j, iter->powers + j, j, j);

    /* bound error after update */
    fmpz_mul2_uiui(iter->pow_error, iter->pow_error,
        iter->max_power, iter->max_power);

    /* readjust precision */
    if (n % 64 == 0 && n > bernoulli_rev_MIN)
    {
        long new_prec, new_max;

        new_prec = global_prec(n);
        new_max = zeta_terms(n, new_prec);

        if (new_prec < iter->prec && new_max <= iter->max_power)
        {
            /* change precision of the powers */
            for (j = 3; j <= new_max; j += 2)
                fmpz_tdiv_q_2exp(iter->powers + j, iter->powers + j,
                    iter->prec - new_prec);

            /* the error also changes precision */
            fmpz_cdiv_q_2exp(iter->pow_error, iter->pow_error,
                iter->prec - new_prec);
            /* contribution of rounding error when changing the precision
               of the powers */
            fmpz_add_ui(iter->pow_error, iter->pow_error, 1);

            /* speed improvement (could be skipped with better multiplication) */
            fmprb_set_round(iter->two_pi_squared, iter->two_pi_squared, new_prec);

            iter->max_power = new_max;
            iter->prec = new_prec;
        }
    }

    iter->n -= 2;

    fmpz_clear(sum);
    fmpr_clear(err);
    fmprb_clear(z);
    fmprb_clear(h);
}

