/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

slong _arb_compute_bs_exponents(slong * tab, slong n);
slong _arb_get_exp_pos(const slong * tab, slong step);

static void
bsplit(fmpz_t T, fmpz_t Q, flint_bitcnt_t * Qexp,
    const slong * xexp,
    const fmpz * xpow, flint_bitcnt_t r, slong a, slong b)
{
    int cc;

    if (b - a == 1)
    {
        count_trailing_zeros(cc, (2 * a + 2));
        fmpz_neg_ui(Q, (2 * a + 2) >> cc);
        fmpz_mul_ui(Q, Q, 2 * a + 3);
        *Qexp = 2 * r + cc;

        fmpz_set(T, xpow);
    }
    else if (b - a == 2)
    {
        fmpz_mul2_uiui(T, xpow, (2 * a + 4), (2 * a + 5));
        fmpz_mul_2exp(T, T, 2 * r);
        fmpz_neg(T, T);
        fmpz_add(T, T, xpow + 1);

        count_trailing_zeros(cc, (2 * a + 4));
        fmpz_neg_ui(Q, (2 * a + 4) >> cc);
        fmpz_mul_ui(Q, Q, 2 * a + 5);
        *Qexp = 2 * r + cc;

        count_trailing_zeros(cc, (2 * a + 2));
        fmpz_mul2_uiui(Q, Q, (2 * a + 2) >> cc, (2 * a + 3));
        fmpz_neg(Q, Q);
        *Qexp += 2 * r + cc;
    }
    else
    {
        slong step, m, i;
        flint_bitcnt_t Q2exp[1];
        fmpz_t Q2, T2;

        step = (b - a) / 2;
        m = a + step;

        fmpz_init(Q2);
        fmpz_init(T2);

        bsplit(T,  Q,  Qexp,  xexp, xpow, r, a, m);
        bsplit(T2, Q2, Q2exp, xexp, xpow, r, m, b);

        fmpz_mul(T, T, Q2);
        fmpz_mul_2exp(T, T, *Q2exp);

        /* find x^step in table */
        i = _arb_get_exp_pos(xexp, step);
        fmpz_addmul(T, xpow + i, T2);  
        fmpz_clear(T2);

        fmpz_mul(Q, Q, Q2);
        *Qexp = *Qexp + *Q2exp;
        fmpz_clear(Q2);
    }
}

/* todo: also allow computing cos, using the same table... */
void
_arb_sin_sum_bs_powtab(fmpz_t T, fmpz_t Q, flint_bitcnt_t * Qexp,
    const fmpz_t x, flint_bitcnt_t r, slong N)
{
    slong * xexp;
    slong length, i;
    fmpz * xpow;

    /* compute the powers of x^2 that will appear (at least x^2) */
    xexp = flint_calloc(2 * FLINT_BITS, sizeof(slong));
    length = _arb_compute_bs_exponents(xexp, N);

    xpow = _fmpz_vec_init(length);
    fmpz_mul(xpow, x, x);

    /* build x^i table */
    for (i = 1; i < length; i++)
    {
        if (xexp[i] == 2 * xexp[i-1])
        {
            fmpz_mul(xpow + i, xpow + i - 1, xpow + i - 1);
        }
        else if (xexp[i] == 2 * xexp[i-2]) /* prefer squaring if possible */
        {
            fmpz_mul(xpow + i, xpow + i - 2, xpow + i - 2);
        }
        else if (xexp[i] == 2 * xexp[i-1] + 1)
        {
            fmpz_mul(xpow + i, xpow + i - 1, xpow + i - 1);
            fmpz_mul(xpow + i, xpow + i, xpow);
        }
        else if (xexp[i] == 2 * xexp[i-2] + 1)
        {
            fmpz_mul(xpow + i, xpow + i - 2, xpow + i - 2);
            fmpz_mul(xpow + i, xpow + i, xpow);
        }
        else
        {
            flint_printf("power table has the wrong structure!\n");
            flint_abort();
        }
    }

    bsplit(T, Q, Qexp, xexp, xpow, r, 0, N);
    _fmpz_vec_clear(xpow, length);
    flint_free(xexp);
}

/*
Determine N such that the error is bounded by 2^-prec when summing the
Taylor series of sin(x) up to term x^(2N+1) inclusive. We choose an N with
many trailing zeros to improve efficiency of the binary splitting.
*/
static slong
bs_num_terms(slong mag, slong prec)
{
    slong N;

    N = _arb_exp_taylor_bound(mag, prec);
    N = N / 2 - 1;
    N = FLINT_MAX(N, 1);

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
arb_sin_cos_fmpz_div_2exp_bsplit(arb_t wsin, arb_t wcos, const fmpz_t x, flint_bitcnt_t r, slong prec)
{
    fmpz_t T, Q;
    slong N, xmag;
    flint_bitcnt_t Qexp[1];

    /* slightly reduce memory usage at very high precision */
    arb_zero(wsin);
    arb_zero(wcos);

    fmpz_init(T);
    fmpz_init(Q);

    if (r > prec)
        flint_abort();

    /* Binary splitting (+1 fixed-point ulp truncation error). */
    xmag = fmpz_bits(x) - r;
    N = bs_num_terms(xmag, prec);
   _arb_sin_sum_bs_powtab(T, Q, Qexp, x, r, N);

    /* we still need to multiply and add x/2^r to get sine */
    fmpz_mul(T, T, x);
    Qexp[0] += r;

    /* T = T / Q  (+1 fixed-point ulp error). */
    if (Qexp[0] >= prec)
    {
        fmpz_tdiv_q_2exp(T, T, Qexp[0] - prec);
        fmpz_tdiv_q(T, T, Q);
    }
    else
    {
        fmpz_mul_2exp(T, T, prec - Qexp[0]);
        fmpz_tdiv_q(T, T, Q);
    }

    fmpz_mul_2exp(Q, x, prec - r);
    fmpz_add(T, T, Q);

    /* T = sin(u) with at most 2 fixed-point ulp error. */
    arf_set_fmpz(arb_midref(wsin), T);
    arf_mul_2exp_si(arb_midref(wsin), arb_midref(wsin), -prec);
    mag_set_ui_2exp_si(arb_radref(wsin), 2, -prec);

    /* compute cos from sin */
    arb_mul(wcos, wsin, wsin, prec);
    arb_sub_ui(wcos, wcos, 1, prec);
    arb_neg(wcos, wcos);
    arb_sqrt(wcos, wcos, prec);

    fmpz_clear(T);
    fmpz_clear(Q);
}

void
arb_sin_cos_arf_bb(arb_t zsin, arb_t zcos, const arf_t x, slong prec)
{
    slong k, iter, bits, r, xmag, q, wp;
    slong argred_bits, start_bits;
    int inexact, negative;
    fmpz_t t, u;
    arb_t wcos, wsin, tmp1;

    if (zsin == NULL)
    {
        arb_init(tmp1);
        arb_sin_cos_arf_bb(tmp1, zcos, x, prec);
        arb_clear(tmp1);
        return;
    }

    if (zcos == NULL)
    {
        arb_init(tmp1);
        arb_sin_cos_arf_bb(zsin, tmp1, x, prec);
        arb_clear(tmp1);
        return;
    }

    if (arf_is_zero(x))
    {
        arb_zero(zsin);
        arb_one(zcos);
        return;
    }

    xmag = arf_abs_bound_lt_2exp_si(x);
    negative = arf_sgn(x) < 0;

    /* We assume that this function only gets called with something
       reasonable as input (huge/tiny input will be handled by
       the main sin/cos wrapper). */
    if (arf_is_special(x) || arf_cmpabs_d(x, 3.15) > 0 || xmag < -2 * prec - 100)
    {
        flint_printf("arb_sin_cos_arf_bb: unexpectedly large/small input\n");
        flint_abort();
    }

    argred_bits = 24;
    start_bits = argred_bits * 3;

    q = FLINT_MAX(0, xmag + argred_bits);
    if (q <= 2)
        q = 0;

    wp = prec + 10 + 2 * (q - xmag) + 2 * FLINT_BIT_COUNT(prec);

    fmpz_init(t);
    fmpz_init(u);
    arb_init(wcos);
    arb_init(wsin);
    arb_init(tmp1);

    /* Convert x/2^q to a fixed-point number. */
    inexact = arf_get_fmpz_fixed_si(t, x, -wp + q);
    fmpz_abs(t, t);

    /* Aliasing of z and x is safe now that only use t. */
    /* Start with z = 1. */
    arb_one(zcos);
    arb_zero(zsin);

    /* Bit-burst loop. */
    for (iter = 0, bits = start_bits; !fmpz_is_zero(t); iter++, bits *= 3)
    {
        /* Extract bits. */
        r = FLINT_MIN(bits, wp);
        fmpz_tdiv_q_2exp(u, t, wp - r);

        arb_sin_cos_fmpz_div_2exp_bsplit(wsin, wcos, u, r, wp);

        /* Remove used bits. */
        fmpz_mul_2exp(u, u, wp - r);
        fmpz_sub(t, t, u);

        /* zsin, zcos = zsin wcos + zcos wsin, zcos wcos - zsin wsin */
        /* using karatsuba */
        arb_add(tmp1, zsin, zcos, wp);
        arb_mul(zcos, zcos, wcos, wp);
        arb_add(wcos, wcos, wsin, wp);
        arb_mul(wsin, wsin, zsin, wp);
        arb_mul(tmp1, tmp1, wcos, wp);
        arb_sub(zsin, tmp1, wsin, wp);
        arb_sub(zsin, zsin, zcos, wp);
        arb_sub(zcos, zcos, wsin, wp);
        arb_zero(tmp1);  /* slightly reduce memory usage */
    }

    /* Initial fixed-point truncation error. */
    if (inexact)
    {
        arb_add_error_2exp_si(zcos, -wp);
        arb_add_error_2exp_si(zsin, -wp);
    }

    if (q != 0)
    {
        /* cos(x) = 2 cos(x/2)^2 - 1 */
        for (k = 0; k < q; k++)
        {
            arb_mul(zcos, zcos, zcos, wp);
            arb_mul_2exp_si(zcos, zcos, 1);
            arb_sub_ui(zcos, zcos, 1, wp);
        }

        arb_mul(tmp1, zcos, zcos, wp);
        arb_sub_ui(tmp1, tmp1, 1, wp);
        arb_neg(tmp1, tmp1);
        arb_sqrt(zsin, tmp1, wp);
    }

    if (negative)
        arb_neg(zsin, zsin);

    arb_set_round(zsin, zsin, prec);
    arb_set_round(zcos, zcos, prec);

    fmpz_clear(t);
    fmpz_clear(u);
    arb_clear(wcos);
    arb_clear(wsin);
    arb_clear(tmp1);
}

