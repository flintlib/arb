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

#include "elefun.h"

/*
Determine N such that the error is bounded by 2^-prec when summing the
Taylor series of exp(x) up to term x^N inclusive. We choose an N with
many trailing zeros to improve efficiency of the binary splitting.
*/
static long
bs_num_terms(long mag, long prec)
{
    long N;

    N = elefun_exp_taylor_bound(mag, prec);
    /* Convert from N exclusive to N inclusive. */
    N--;

    if (N > 10000)
        while (N % 128 != 0)
            N++;

    if (N > 1000)
        while (N % 16 != 0)
            N++;

    if (N > 100)
        while (N % 2 != 0)
            N++;

    return N;
}

void
elefun_exp_fmpr_bb(fmprb_t z, const fmpr_t x, long prec, int m1)
{
    long k, iter, bits, r, mag, q, wp, N;
    long argred_bits, start_bits;
    mp_bitcnt_t Qexp[1];
    int inexact;
    fmpz_t t, u, T, Q;

    if (fmpr_is_zero(x))
    {
        fmprb_one(z);
        return;
    }

    mag = fmpr_abs_bound_lt_2exp_si(x);

    /* We assume that this function only gets called with something
       reasonable as input (huge/tiny input will be handled by
       the main exp wrapper). */
    if (mag > 200 || mag < -2 * prec - 100)
    {
        printf("elefun_exp_fmpr_bb: unexpectedly large/small input\n");
        abort();
    }

    if (prec < 100000000)
    {
        argred_bits = 16;
        start_bits = 32;
    }
    else
    {
        argred_bits = 32;
        start_bits = 64;
    }

    /* Argument reduction: exp(x) -> exp(x/2^q). This improves efficiency
       of the first iteration in the bit-burst algorithm. */
    q = FLINT_MAX(0, mag + argred_bits);

    /* Determine working precision. */
    wp = prec + 10 + 2 * q + 2 * FLINT_BIT_COUNT(prec);
    if (m1 && mag < 0)
        wp += (-mag);

    fmpz_init(t);
    fmpz_init(u);
    fmpz_init(Q);
    fmpz_init(T);

    /* Convert x/2^q to a fixed-point number. */
    inexact = fmpr_get_fmpz_fixed_si(t, x, -wp + q);

    /* Aliasing of z and x is safe now that only use t. */
    /* Start with z = 1. */
    fmprb_one(z);

    /* Bit-burst loop. */
    for (iter = 0, bits = start_bits; !fmpz_is_zero(t);
        iter++, bits *= 2)
    {
        /* Extract bits. */
        r = FLINT_MIN(bits, wp);
        fmpz_tdiv_q_2exp(u, t, wp - r);

        /* Binary splitting (+1 fixed-point ulp truncation error). */
        mag = fmpz_bits(u) - r;
        N = bs_num_terms(mag, wp);
        elefun_exp_sum_bs_powtab(T, Q, Qexp, u, r, N);

        /* T = T / Q  (+1 fixed-point ulp error). */
        if (*Qexp >= wp)
        {
            fmpz_tdiv_q_2exp(T, T, *Qexp - wp);
            fmpz_tdiv_q(T, T, Q);
        }
        else
        {
            fmpz_mul_2exp(T, T, wp - *Qexp);
            fmpz_tdiv_q(T, T, Q);
        }

        /* T = 1 + T */
        fmpz_one(Q);
        fmpz_mul_2exp(Q, Q, wp);
        fmpz_add(T, T, Q);

        /* Now T = exp(u) with at most 2 fixed-point ulp error. */
        /* Set z = z * T. */
        {
            fmprb_t w;
            fmprb_init(w);
            fmpr_set_fmpz(fmprb_midref(w), T);
            fmpr_mul_2exp_si(fmprb_midref(w), fmprb_midref(w), -wp);
            fmpr_set_si_2exp_si(fmprb_radref(w), 2, -wp);
            fmprb_mul(z, z, w, wp);
            fmprb_clear(w);
        }

        /* Remove used bits. */
        fmpz_mul_2exp(u, u, wp - r);
        fmpz_sub(t, t, u);
    }

    /* We have exp(x + eps) - exp(x) < 2*eps (by assumption that the argument
       reduction is large enough). */
    if (inexact)
        fmprb_add_error_2exp_si(z, -wp + 1);

    fmpz_clear(t);
    fmpz_clear(u);
    fmpz_clear(Q);
    fmpz_clear(T);

    /* exp(x) = exp(x/2^q)^(2^q) */
    for (k = 0; k < q; k++)
        fmprb_mul(z, z, z, wp);

    if (m1)
        fmprb_sub_ui(z, z, 1, wp);

    fmprb_set_round(z, z, prec);
}

