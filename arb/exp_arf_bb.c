/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

/*
Determine N such that the error is bounded by 2^-prec when summing the
Taylor series of exp(x) up to term x^N inclusive. We choose an N with
many trailing zeros to improve efficiency of the binary splitting.
*/
static slong
bs_num_terms(slong mag, slong prec)
{
    slong N;

    N = _arb_exp_taylor_bound(mag, prec);
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
arb_exp_arf_bb(arb_t z, const arf_t x, slong prec, int minus_one)
{
    slong k, iter, bits, r, mag, q, wp, N;
    slong argred_bits, start_bits;
    flint_bitcnt_t Qexp[1];
    int inexact;
    fmpz_t t, u, T, Q;
    arb_t w;

    if (arf_is_zero(x))
    {
        if (minus_one)
            arb_zero(z);
        else
            arb_one(z);
        return;
    }

    if (arf_is_special(x))
    {
        flint_abort();
    }

    mag = arf_abs_bound_lt_2exp_si(x);

    /* We assume that this function only gets called with something
       reasonable as input (huge/tiny input will be handled by
       the main exp wrapper). */
    if (mag > 200 || mag < -2 * prec - 100)
    {
        flint_printf("arb_exp_arf_bb: unexpectedly large/small input\n");
        flint_abort();
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
    if (minus_one && mag < 0)
        wp += (-mag);

    fmpz_init(t);
    fmpz_init(u);
    fmpz_init(Q);
    fmpz_init(T);
    arb_init(w);

    /* Convert x/2^q to a fixed-point number. */
    inexact = arf_get_fmpz_fixed_si(t, x, -wp + q);

    /* Aliasing of z and x is safe now that only use t. */
    /* Start with z = 1. */
    arb_one(z);

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

       _arb_exp_sum_bs_powtab(T, Q, Qexp, u, r, N);

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
        arf_set_fmpz(arb_midref(w), T);
        arf_mul_2exp_si(arb_midref(w), arb_midref(w), -wp);
        mag_set_ui_2exp_si(arb_radref(w), 2, -wp);
        arb_mul(z, z, w, wp);

        /* Remove used bits. */
        fmpz_mul_2exp(u, u, wp - r);
        fmpz_sub(t, t, u);
    }

    /* We have exp(x + eps) - exp(x) < 2*eps (by assumption that the argument
       reduction is large enough). */
    if (inexact)
        arb_add_error_2exp_si(z, -wp + 1);

    fmpz_clear(t);
    fmpz_clear(u);
    fmpz_clear(Q);
    fmpz_clear(T);
    arb_clear(w);

    /* exp(x) = exp(x/2^q)^(2^q) */
    for (k = 0; k < q; k++)
        arb_mul(z, z, z, wp);

    if (minus_one)
        arb_sub_ui(z, z, 1, wp);

    arb_set_round(z, z, prec);
}

