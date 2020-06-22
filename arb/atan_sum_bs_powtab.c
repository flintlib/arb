/*
    Copyright (C) 2014 Fredrik Johansson

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
    if (b - a == 1)
    {
        fmpz_set(T, xpow);

        if (a % 2 == 0)
            fmpz_neg_ui(Q, 2 * a + 3);
        else
            fmpz_set_ui(Q, 2 * a + 3);

        *Qexp = 2 * r;
    }
    else if (b - a == 2)
    {
        fmpz_mul_ui(T, xpow, 2 * a + 5);
        fmpz_mul_2exp(T, T, 2 * r);
        fmpz_submul_ui(T, xpow + 1, 2 * a + 3);

        if (a % 2 == 1)
            fmpz_neg(T, T);

        fmpz_neg_ui(Q, 2 * a + 3);
        fmpz_mul_ui(Q, Q, 2 * a + 5);
        *Qexp = 4 * r;
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
        fmpz_mul(T2, T2, Q);
        fmpz_addmul(T, xpow + i, T2);
        fmpz_clear(T2);

        fmpz_mul(Q, Q, Q2);
        *Qexp = *Qexp + *Q2exp;
        fmpz_clear(Q2);
    }
}

void
_arb_atan_sum_bs_powtab(fmpz_t T, fmpz_t Q, flint_bitcnt_t * Qexp,
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

