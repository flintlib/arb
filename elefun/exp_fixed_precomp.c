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

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include "elefun.h"

/* todo: dynamic alloc */
__thread fmpz_t ln2_cache;
__thread fmpz exp_cache[EXP_CACHE_LEVELS][EXP_CACHE_NUM];
__thread long exp_cache_exp[EXP_CACHE_LEVELS][EXP_CACHE_NUM];
__thread int exp_cache_init = 0;

/* max 1 ulp error (TODO: check this) */
void
compute_exp_cache()
{
    long i, j, wp;
    fmpr_t x, y;

    wp = EXP_CACHE_PREC + 30;

    fmpr_init(x);
    fmpr_init(y);

    for (i = 0; i < EXP_CACHE_LEVELS; i++)
    {
        fmpr_set_ui_2exp_si(x, 1, -(i+1) * EXP_CACHE_BITS);
        fmpr_exp(x, x, wp, FMPR_RND_DOWN);
        fmpr_one(y);

        for (j = 0; j < EXP_CACHE_NUM; j++)
        {
            fmpz_init(exp_cache[i] + j);
            fmpr_get_fmpz_fixed_si(exp_cache[i] + j, y, -EXP_CACHE_PREC);
            fmpr_mul(y, y, x, wp, FMPR_RND_DOWN);
        }
    }

    fmpz_init(ln2_cache);
    fmpr_set_ui(x, 2);
    fmpr_log(x, x, wp, FMPR_RND_DOWN);

    fmpr_get_fmpz_fixed_si(ln2_cache, x, -EXP_CACHE_PREC);

    fmpr_clear(x);
    fmpr_clear(y);

    exp_cache_init = 1;
}

void
elefun_exp_fixed_precomp(fmpz_t y, fmpz_t yerr, fmpz_t exponent,
    const fmpz_t x, const fmpz_t xerr, long prec)
{
    long i, n;
    ulong r, red;
    fmpz_t t, u, uerr;

    fmpz_init(t);
    fmpz_init(u);
    fmpz_init(uerr);

    if (!exp_cache_init)
        compute_exp_cache();

    /* x = u + exponent * ln2 */
    fmpz_tdiv_q_2exp(t, ln2_cache, EXP_CACHE_PREC - prec);
    fmpz_fdiv_qr(exponent, u, x, t);

    /* ln2 has 1 ulp error -> u gains exponent * ulp error */
    if (fmpz_sgn(exponent) >= 0)
        fmpz_add(uerr, xerr, exponent);
    else
        fmpz_sub(uerr, xerr, exponent);

    /* extract top bits of u */
    if (prec < EXP_CACHE_LEVELS * EXP_CACHE_BITS)
    {
        fmpz_mul_2exp(t, u, EXP_CACHE_LEVELS * EXP_CACHE_BITS - prec);
        red = fmpz_get_ui(t);
        fmpz_zero(u);
    }
    else
    {
        fmpz_tdiv_q_2exp(t, u, prec - EXP_CACHE_LEVELS * EXP_CACHE_BITS);
        red = fmpz_get_ui(t);
        fmpz_mul_2exp(t, t, prec - EXP_CACHE_LEVELS * EXP_CACHE_BITS);
        fmpz_sub(u, u, t);
    }

    /* Taylor series for remaining bits */
    n = elefun_exp_taylor_bound(-EXP_CACHE_REDUCTION, prec);
    elefun_exp_fixed_taylor_horner_precomp(y, yerr, u, n, prec);

    /* 1 ulp error for truncation of Taylor series */
    fmpz_add_ui(yerr, yerr, 1);

    /* cache for high bits */
    for (i = EXP_CACHE_LEVELS - 1; i >= 0; i--)
    {
        r = (red >> ((EXP_CACHE_LEVELS - i - 1) * EXP_CACHE_BITS));
        r &= ((1UL << (EXP_CACHE_BITS)) - 1);

        fmpz_tdiv_q_2exp(u, exp_cache[i] + r, EXP_CACHE_PREC - prec);
        fmpz_mul_tdiv_q_2exp(y, y, u, prec);
    }

    if (EXP_CACHE_LEVELS == 2)
    {
        fmpz_mul_ui(yerr, yerr, 7);
        fmpz_add_ui(yerr, yerr, 9);
    }
    else
        abort();

    /* propagated error <= (err(u) * exp(uerr)) * exp(x) */

    /* set uerr to a bound for err(u) * exp(err(u)) */
    if (fmpz_bits(uerr) <= prec)   /* for a < 1.79329, a*exp(a) < a + a^2 + a^3 */
    {
        fmpz_mul(u, uerr, uerr);  /* uerr^2 */
        fmpz_cdiv_q_2exp(u, u, prec);

        fmpz_mul(t, u, uerr);     /* uerr^3 */
        fmpz_cdiv_q_2exp(t, t, prec);

        fmpz_add(uerr, uerr, t);
        fmpz_add(uerr, uerr, u);
    }
    else
    {
        /* crude bound: a*exp(a) <= 4^ceil(a) */
        fmpz_cdiv_q_2exp(u, uerr, prec);

        if (fmpz_bits(u) > 20)
            abort();

        fmpz_one(uerr);
        fmpz_mul_2exp(uerr, uerr, prec + 2 * (*u));
    }

    fmpz_add(t, y, yerr); /* upper bound for exp(u) */

    fmpz_mul(uerr, uerr, t); /* multiply by (err(u) * exp(err(u))) */
    fmpz_cdiv_q_2exp(uerr, uerr, prec);

    fmpz_add(yerr, yerr, uerr);  /* add to final error for y = exp(x) */

    fmpz_clear(t);
    fmpz_clear(u);
    fmpz_clear(uerr);
}

