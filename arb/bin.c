/*
    Copyright (C) 2013 Fredrik Johansson
    Copyright (C) 2021, 2022 Albin Ahlbäck

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "flint/fmpz.h"

#define LOG2_BIN(n, k)                      \
    (-1.3257480647361592                    \
     + ((n) + .5) * log2(n)                 \
     - ((k) + .5) * log2(k)                 \
     - ((n) - (k) + .5) * log2((n) - (k)))  \

static void
_arb_bin_ui(arb_t x, const arb_t n, ulong k, slong prec)
{
    arb_t t, u;

    arb_init(t);
    arb_init(u);

    arb_sub_ui(t, n, k - 1, prec);
    arb_rising_ui(t, t, k, prec);
    arb_fac_ui(u, k, prec);
    arb_div(x, t, u, prec);

    arb_clear(t);
    arb_clear(u);
}

void
arb_bin_ui(arb_t x, const arb_t n, ulong k, slong prec)
{
    if (k == 0)
        arb_one(x);
    else if (k == 1)
        arb_set_round(x, n, prec);
    else if (arb_is_int(n)
            && arb_is_nonnegative(n)
            && ARF_EXP(arb_midref(n)) <= FLINT_BITS) /* fits in ulong */
    {
        fmpz_t tmp;
        fmpz_init(tmp);
        arb_get_unique_fmpz(tmp, n);
        arb_bin_uiui(x, fmpz_get_ui(tmp), k, prec);
    }
    else
        _arb_bin_ui(x, n, k, prec);
}

void
arb_bin_uiui(arb_t x, ulong n, ulong k, slong prec)
{
    k = FLINT_MIN(k, n - k);

#if LONG_MAX == WORD_MAX
    if (n <= (UWORD(1) << 12)
            || k <= (UWORD(1) << 7)
#if FLINT_BITS == 64
            || (n <= (UWORD(1) << 32) && k <= (UWORD(1) << 8))
#else
            || (n <= UWORD_MAX && k <= (UWORD(1) << 8))
#endif
            || (n <= (UWORD(1) << 22) && k <= (UWORD(1) << 9))
            || (n <= (UWORD(1) << 16) && k <= (UWORD(1) << 10))
            || ((double) prec) < LOG2_BIN(n, k))
    {
        mpz_t mres;
        mpz_init(mres);
        flint_mpz_bin_uiui(mres, n, k);
        arb_set_round_mpz(x, mres, prec);
        mpz_clear(mres);
    }
#else
    if (n <= (UWORD(1) << 12)
            || k <= (UWORD(1) << 7)
            || (n < UWORD_MAX && k <= (UWORD(1) << 8))
            || (n <= (UWORD(1) << 22) && k <= (UWORD(1) << 9))
            || (n <= (UWORD(1) << 16) && k <= (UWORD(1) << 10)))
    {
        mpz_t mres;
        mpz_init(mres);
        flint_mpz_bin_uiui(mres, n, k);
        arb_set_round_mpz(x, mres, prec);
        mpz_clear(mres);
    }
    else if (k < (UWORD(1) << 32) && ((double) prec) < LOG2_BIN(n, k))
    {
        mpz_t mres;
        __mpz_struct mn = {1, 1, NULL};
        mpz_init(mres);
        mn._mp_d = &n;
        mpz_bin_ui(mres, &mn, k);
        arb_set_round_mpz(x, mres, prec);
        mpz_clear(mres);
    }
#endif
    else
    {
        arf_struct atmp;
        mag_struct mtmp = {0, 0};
        arb_struct tmp = {atmp, mtmp};
        arf_init_set_ui(&atmp, n);
        arb_bin_ui(x, &tmp, k, prec);
    }
}

#undef LOG2_BIN
